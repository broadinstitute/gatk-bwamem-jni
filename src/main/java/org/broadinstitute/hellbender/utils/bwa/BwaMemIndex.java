package org.broadinstitute.hellbender.utils.bwa;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * BwaMemIndex manages the mapping of a bwa index image file into (non-Java) memory.
 * It's typically a huge chunk of memory, so you need to manage it as a precious resource.
 *
 * Usage pattern is:
 *   Create a BwaMemIndex for some reference.
 *   Create Aligners as needed to do some aligning -- they're pretty lightweight and thread safe.
 *     (But you may need to manage memory by controlling the number of reads that you align in one chunk.)
 *   Close the BwaMemIndex when you're done aligning.
 *
 * This class doesn't know anything about Spark.  You can use it in a distributed setting if you distribute the index
 * file to each node using the Yarn --files mechanism.  You might find it convenient to manage a singleton instance
 * of a BwaMemIndex on each Java VM when you're running distributed:  check out BwaMemIndexSingleton in GATK.
 *
 * Alternatively, you could use this class directly to run bwa multi-threaded on a single machine,
 * if that's what you want to do.
 */
public final class BwaMemIndex implements AutoCloseable {
    private final String indexImageFile; // stash this for error messages
    private volatile long indexAddress; // address where the index was memory-mapped (for use by C code)
    private final AtomicInteger refCount; // keep track of how many threads are actively aligning
    private final List<String> refContigNames; // the reference dictionary from the index
    private static volatile boolean nativeLibLoaded = false; // whether we've loaded the native library or not

    /**
     * This method creates the new, single-file bwa index image from the 5 files that the "bwa index" command produces.
     * It throws an IllegalStateException if it is unable to create the index image.
     */
    public static void createIndexImage( final String referenceName, final String imageName ) {
        if ( referenceName == null ) {
            throw new IllegalStateException("referenceName cannot be null.");
        }
        if ( imageName == null ) {
            throw new IllegalStateException("imageName cannot be null.");
        }
        loadNativeLibrary();
        assertNonEmptyReadable(referenceName+".amb");
        assertNonEmptyReadable(referenceName+".ann");
        assertNonEmptyReadable(referenceName+".bwt");
        assertNonEmptyReadable(referenceName+".pac");
        assertNonEmptyReadable(referenceName+".sa");
        String errMsg = createIndexImageFile(referenceName, imageName);
        if ( errMsg != null ) {
            throw new IllegalStateException(errMsg);
        }
    }

    /** create an index from an image file.  use BwaMemIndex.createIndexImage to create one, as necessary. */
    public BwaMemIndex( final String indexImageFile ) {
        this(indexImageFile, false, false);
    }

    public BwaMemIndex( final String indexImageFile, final boolean ignoreVersion, final boolean compareCRC ) {
        if ( indexImageFile == null ) {
            throw new IllegalStateException("indexImageFile cannot be null.");
        }
        this.indexImageFile = indexImageFile;
        loadNativeLibrary();
        assertNonEmptyReadable(indexImageFile);
        refCount = new AtomicInteger();
        final String addrOrErrorMessage = openIndex(indexImageFile, ignoreVersion, compareCRC);
        try {
            indexAddress = Long.parseLong(addrOrErrorMessage);
        } catch ( NumberFormatException nfe ) {
            throw new IllegalStateException(addrOrErrorMessage);
        }

        final ByteBuffer refContigNamesBuf = getRefContigNames(indexAddress);
        if ( refContigNamesBuf == null ) {
            throw new IllegalStateException("Unable to retrieve reference contig names from bwa-mem index "+indexImageFile);
        }
        refContigNamesBuf.order(ByteOrder.nativeOrder()).position(0).limit(refContigNamesBuf.capacity());
        final int nRefContigNames = refContigNamesBuf.getInt();
        refContigNames = new ArrayList<>(nRefContigNames);
        for ( int idx = 0; idx < nRefContigNames; ++idx ) {
            int nameLen = refContigNamesBuf.getInt();
            byte[] nameBytes = new byte[nameLen];
            refContigNamesBuf.get(nameBytes);
            refContigNames.add(new String(nameBytes));
        }
        destroyByteBuffer(refContigNamesBuf);
    }

    /** true if index has not been closed */
    public boolean isOpen() { return indexAddress != 0L; }

    /** there's someone using the index -- don't allow it to be closed */
    public long refIndex() {
        refCount.incrementAndGet();
        if ( indexAddress == 0L ) {
            throw new IllegalStateException("Index image " + indexImageFile + " has been closed");
        }
        return indexAddress;
    }

    /** done using the index -- if ref count has fallen to 0, a call to close can be expected to succeed */
    public void deRefIndex() { refCount.decrementAndGet(); }

    /** close the index and release the (non-Java) memory that's been allocated */
    @Override
    public void close() {
        long addr = indexAddress;
        if ( addr != 0L ) {
            synchronized (BwaMemIndex.class) {
                addr = indexAddress;
                if ( addr != 0L ) {
                    if ( refCount.intValue() != 0 ) {
                        throw new IllegalStateException("Index image "+indexImageFile+" can't be closed:  it's in use.");
                    }
                    indexAddress = 0L;
                    final String errMsg = destroyIndex(indexAddress);
                    if ( errMsg != null ) {
                        throw new IllegalStateException(indexImageFile+" "+errMsg);
                    }
                }
            }
        }
    }

    /** retrieve list of contig names in the reference dictionary */
    public List<String> getReferenceContigNames() {
        return refContigNames;
    }

    /** returns github GUID for the version of bwa that has been compiled */
    public static String getBWAVersion() {
        loadNativeLibrary();
        return getVersion();
    }

    ByteBuffer doAlignment( final ByteBuffer seqs, final ByteBuffer opts ) {
        final ByteBuffer alignments = createAlignments(seqs, indexAddress, opts);
        if ( alignments == null ) {
            throw new IllegalStateException("Unable to get alignments from bwa-mem index "+indexImageFile+": We don't know why.");
        }
        return alignments;
    }

    private static void assertNonEmptyReadable( final String fileName ) {
        if ( !nonEmptyReadableFile(fileName) )
            throw new IllegalArgumentException("Missing bwa index file: "+fileName);
    }

    private static boolean nonEmptyReadableFile( final String fileName ) {
        try ( final FileInputStream is = new FileInputStream(fileName) ) {
            return is.read() != -1;
        } catch ( final IOException ioe ) {
            return false;
        }
    }

    private static void loadNativeLibrary() {
        if ( !nativeLibLoaded ) {
            synchronized(BwaMemIndex.class) {
                if ( !nativeLibLoaded ) {
                    final String libNameOverride = System.getProperty("LIBBWA_PATH");
                    if ( libNameOverride != null ) {
                        System.load(libNameOverride);
                    }
                    else {
                        final String osName = System.getProperty("os.name", "unknown").toUpperCase();
                        final String osArch = System.getProperty("os.arch");
                        final String libName;
                        if ( !"x86_64".equals(osArch) && !"amd64".equals(osArch) ) {
                            throw new IllegalStateException(
                                    "We have pre-built fermi-lite binaries only for x86_64 and amd64.  "+
                                            "Your os.arch is "+osArch+"."+
                                            "Set property LIBBWA_PATH to point to a native library for your architecture.");
                        }
                        if ( osName.startsWith("MAC") ) libName = "/libbwa.Darwin.dylib";
                        else if ( osName.startsWith("LINUX") ) libName = "/libbwa.Linux.so";
                        else {
                            throw new IllegalStateException(
                                    "We have pre-built fermi-lite binaries only for Linux and Mac.  "+
                                            "Your os.name is "+osName+"."+
                                            "Set property LIBBWA_PATH to point to a native library for your operating system.");
                        }
                        try ( final InputStream is = BwaMemIndex.class.getResourceAsStream(libName) ) {
                            if ( is == null ) {
                                throw new IllegalStateException("Can't find resource "+libName);
                            }
                            final File tmpFile = File.createTempFile("libbwa.",".jnilib");
                            tmpFile.deleteOnExit();
                            Files.copy(is, tmpFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                            System.load(tmpFile.getPath());
                        }
                        catch (IOException ioe ) {
                            throw new IllegalStateException("Misconfiguration: Unable to load fermi-lite native library "+libName, ioe);
                        }
                    }
                    nativeLibLoaded = true;
                }
            }
        }
    }

    private static native String createIndexImageFile( String referenceName, String imageName );
    private static native String openIndex( String imageName, boolean ignoreVersion, boolean compareCRC );
    private static native String destroyIndex( long indexAddress );
    static native ByteBuffer createDefaultOptions();
    private static native ByteBuffer getRefContigNames( long indexAddress );
    private static native ByteBuffer createAlignments( ByteBuffer seqs, long indexAddress, ByteBuffer opts );
    static native void destroyByteBuffer( ByteBuffer alignments );
    private static native String getVersion();
}
