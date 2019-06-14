package org.broadinstitute.hellbender.utils.bwa;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * BwaMemIndex manages the mapping of a bwa index image file into (non-Java) memory.
 * It's typically a huge chunk of memory, so you need to manage it as a precious resource.
 * <p>
 * Usage pattern is:
 * Create a BwaMemIndex for some reference.
 * Create Aligners as needed to do some aligning -- they're pretty lightweight and thread safe.
 * (But you may need to manage memory by controlling the number of reads that you align in one chunk.)
 * Close the Aligners when you're done aligning.
 * Close the BwaMemIndex when you're all done aligning.
 * <p>
 * This class doesn't know anything about Spark.  You can use it in a distributed setting if you distribute the index
 * file to each node using the Yarn --files mechanism.  You might find it convenient to manage a singleton instance
 * of a BwaMemIndex on each Java VM when you're running distributed:  check out BwaMemIndexSingleton in GATK.
 * <p>
 * Alternatively, you could use this class directly to run bwa multi-threaded on a single machine,
 * if that's what you want to do.
 */
public final class BwaMemIndex implements AutoCloseable {

    public static final List<String> INDEX_FILE_EXTENSIONS =
            Collections.unmodifiableList(Arrays.asList(".amb", ".ann", ".bwt", ".pac", ".sa"));

    public static final String IMAGE_FILE_EXTENSION = ".img";

    private final String indexImageFile; // stash this for error messages
    private volatile long indexAddress; // address where the index was memory-mapped (for use by C code)
    private final AtomicInteger refCount; // keep track of how many threads are actively aligning
    private final List<String> refContigNames; // the reference dictionary from the index
    private static volatile boolean nativeLibLoaded = false; // whether we've loaded the native library or not

    /**
     * Indexing algorithms supported by Bwa.
     */
    public enum Algorithm {

        /**
         * Chooses the most appropriate algorithm based on characteristic
         * of the reference. For references shorter than 2Gbases it will use {@link #IS}
         * whereas for larger reference it will employ {@link #RB2}.
         */
        AUTO,

        /**
         * The <i>Ropebwt2</i> algorithm.
         *
         * @see "https://arxiv.org/pdf/1406.0426.pdf"
         */
        RB2,

        /**
         * code adapted from original BWT paper -- not recommended
         */
        BWTSW,

        /**
         * Linear-time algorithm for constructing suffix array.
         * It requires 5.37N memory where N is the size of the database.
         * IS is moderately fast, but does not work with database larger than 2GB.
         * IS is the default algorithm due to its simplicity. The current codes for
         * IS algorithm are reimplemented by Yuta Mori.
         *
         * @see "http://bio-bwa.sourceforge.net/bwa.shtml"
         */
        IS;
    }

    /**
     * Create the index image file for a complete set of BWA index files
     *
     * @param indexPrefix the location of the index files.
     * @param imageFile   the location of the new index image file.
     * @throws BwaMemException if there is a problem creating the image file.
     *
     * <p>
     * <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt
     *   process interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     */
    public static void createIndexImageFromIndexFiles( final String indexPrefix, final String imageFile ) {
        if ( indexPrefix == null ) {
            throw new IllegalArgumentException("the index prefix cannot be null");
        }
        if ( imageFile == null ) {
            throw new IllegalArgumentException("the image file name cannot be null");
        }
        for ( final String ext : INDEX_FILE_EXTENSIONS ) {
            final String fileName = indexPrefix + ext;
            if ( !nonEmptyReadableFile(fileName) )
                throw new BwaMemException("can't read bwa index file: " + fileName);
        }

        loadNativeLibrary();
        createIndexImageFile(indexPrefix, imageFile);
    }

    /**
     * <p>
     * <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt
     *   process interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     *
     * @param fasta the location of the fasta reference file.
     * @return location of the generated index image file (fastaFile + ".img")
     * @throws BwaMemException if there is a problem indexing the reference, or creating the image file.
     */
    public static String createIndexImageFromFastaFile( final String fasta ) {
        final String imageFile = fasta + IMAGE_FILE_EXTENSION;
        createIndexImageFromFastaFile(fasta, imageFile);
        return imageFile;
    }

    /**
     * Creates the index image file for a reference fasta file.
     * <p>
     * The index will be created using the default algorithm {@link Algorithm#AUTO}.
     * </p>
     * <p>
     * <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt
     *   process interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     *
     * @param fasta     the location of the targeted reference.
     * @param imageFile the location of the new index image file.
     * @throws BwaMemException if there is a problem indexing the reference, or creating the image file.
     */
    public static void createIndexImageFromFastaFile( final String fasta, final String imageFile ) {
        createIndexImageFromFastaFile(fasta, imageFile, Algorithm.AUTO);
    }

    /**
     * Creates the index image file for a reference fasta file.
     * <p>
     * <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt
     *   process interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     *
     * @param fasta     the location of the targeted reference.
     * @param imageFile the location of the new index image file.
     * @param algorithm the algorithm to use to construct the index (see {@link Algorithm} to see what's available.).
     * @throws BwaMemException if there is a problem indexing the reference, or creating the image file.
     */
    public static void createIndexImageFromFastaFile( final String fasta,
                                                      final String imageFile,
                                                      final Algorithm algorithm ) {
        if ( fasta == null ) {
            throw new IllegalArgumentException("the reference fasta file name cannot be null");
        }
        if ( imageFile == null ) {
            throw new IllegalArgumentException("the image file name cannot be null");
        }
        if ( algorithm == null ) {
            throw new IllegalArgumentException("the input algorithm must not be null");
        }

        final String indexPrefix;
        try {
            indexPrefix = File.createTempFile("bwaIndex", "").getPath();
        } catch ( final IOException ioe ) {
            throw new BwaMemException("unable to create temp file for bwa index file naming", ioe);
        }

        loadNativeLibrary();
        indexReference(fasta, indexPrefix, algorithm.ordinal());
        createIndexImageFile(indexPrefix, imageFile);

        for ( final String ext : INDEX_FILE_EXTENSIONS )
            new File(indexPrefix + ext).delete();
        new File(indexPrefix).delete();
    }

    /**
     * Loads an index from an image file.
     * <p>
     * You can use other methods to create such
     * indexes from fasta reference ({@link #createIndexImageFromFastaFile} or their index files
     * ({@link #createIndexImageFromIndexFiles}).
     * </p>
     * <p>
     * <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt
     *   process interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     *
     * @throws BwaMemException   if some problem occurred when loading the image file.
     */
    public BwaMemIndex( final String indexImageFile ) {
        if ( indexImageFile == null ) {
            throw new IllegalArgumentException("indexImageFile is null");
        }
        if ( !nonEmptyReadableFile(indexImageFile) ) {
            throw new BwaMemException("indexImageFile is missing or empty: " + indexImageFile);
        }
        loadNativeLibrary();
        this.indexImageFile = indexImageFile;
        refCount = new AtomicInteger();
        indexAddress = openIndex(indexImageFile);
        if ( indexAddress == 0L ) {
            throw new BwaMemException("unable to open bwa-mem index: " + indexImageFile);
        }
        ByteBuffer refContigNamesBuf = getRefContigNames(indexAddress);
        if ( refContigNamesBuf == null ) {
            throw new BwaMemException("unable to retrieve reference contig names from bwa-mem index");
        }
        refContigNamesBuf.order(ByteOrder.nativeOrder()).position(0).limit(refContigNamesBuf.capacity());
        int nRefContigNames = refContigNamesBuf.getInt();
        refContigNames = new ArrayList<>(nRefContigNames);
        for ( int idx = 0; idx < nRefContigNames; ++idx ) {
            int nameLen = refContigNamesBuf.getInt();
            byte[] nameBytes = new byte[nameLen];
            refContigNamesBuf.get(nameBytes);
            refContigNames.add(new String(nameBytes));
        }
        destroyByteBuffer(refContigNamesBuf);
    }

    /**
     * true if index has not been closed
     */
    public boolean isOpen() {
        return indexAddress != 0L;
    }

    /**
     * Close the index and release the (non-Java) memory that's been allocated
     *
     * <p>
     * <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt
     *   process interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     */
    @Override
    public void close() {
        long addr = indexAddress;
        if ( addr != 0L ) {
            synchronized ( BwaMemIndex.class ) {
                addr = indexAddress;
                if ( addr != 0L ) {
                    if ( refCount.intValue() != 0 ) {
                        throw new BwaMemException("Index image " + indexImageFile + " can't be closed:  it's in use.");
                    }
                    indexAddress = 0L;
                    destroyIndex(addr);
                }
            }
        }
    }

    /**
     * retrieve list of contig names in the reference dictionary
     */
    public List<String> getReferenceContigNames() {
        return refContigNames;
    }

    /**
     * returns github GUID for the version of bwa that has been compiled
     */
    public static String getBWAVersion() {
        loadNativeLibrary();
        final String version = getVersion();
        if ( version == null ) {
            throw new BwaMemException("unable to retrieve version information");
        }
        return version;
    }

    /**
     * there's someone using the index -- don't allow it to be closed
     */
    void refIndex() {
        refCount.incrementAndGet();
        if ( indexAddress == 0L ) {
            throw new BwaMemException("Index image " + indexImageFile + " has been closed");
        }
    }

    /**
     * done using the index -- if ref count has fallen to 0, a call to close can be expected to succeed
     */
    void deRefIndex() {
        refCount.decrementAndGet();
    }

    ByteBuffer doAlignment( final ByteBuffer seqs, final ByteBuffer opts ) {
        final ByteBuffer alignments = createAlignments(seqs, indexAddress, opts);
        if ( alignments == null ) {
            throw new BwaMemException(
                    "Unable to get alignments from bwa-mem index " + indexImageFile + ": We don't know why.");
        }
        return alignments;
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
            synchronized ( BwaMemIndex.class ) {
                if ( !nativeLibLoaded ) {
                    final String libNameOverride = System.getProperty("LIBBWA_PATH");
                    if ( libNameOverride != null ) {
                        System.load(libNameOverride);
                    } else {
                        final String osName = System.getProperty("os.name", "unknown").toUpperCase();
                        final String osArch = System.getProperty("os.arch");
                        final String libName;
                        if ( !"x86_64".equals(osArch) && !"amd64".equals(osArch) ) {
                            throw new BwaMemException(
                                    "We have pre-built fermi-lite binaries only for x86_64 and amd64.  " +
                                    "Your os.arch is " + osArch + "." +
                                    "Set property LIBBWA_PATH to point to a native library for your architecture.");
                        }
                        if ( osName.startsWith("MAC") ) libName = "/libbwa.Darwin.dylib";
                        else if ( osName.startsWith("LINUX") ) libName = "/libbwa.Linux.so";
                        else {
                            throw new BwaMemException(
                                    "We have pre-built fermi-lite binaries only for Linux and Mac.  " +
                                    "Your os.name is " + osName + "." +
                                    "Set property LIBBWA_PATH to point to a native library for your operating system.");
                        }
                        try ( final InputStream is = BwaMemIndex.class.getResourceAsStream(libName) ) {
                            if ( is == null ) {
                                throw new BwaMemException("Can't find resource " + libName);
                            }
                            final File tmpFile = File.createTempFile("libbwa.", ".jnilib");
                            tmpFile.deleteOnExit();
                            Files.copy(is, tmpFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
                            System.load(tmpFile.getPath());
                        } catch ( IOException ioe ) {
                            throw new BwaMemException(
                                    "Misconfiguration: Unable to load fermi-lite native library " + libName, ioe);
                        }
                    }
                    nativeLibLoaded = true;
                }
            }
        }
    }

    private static native void indexReference( String referenceName, String indexPrefix, int algorithm );
    private static native void createIndexImageFile( String indexPrefix, String imageName );
    private static native long openIndex( String indexImageFile );
    private static native void destroyIndex( long indexAddress );
    static native ByteBuffer createDefaultOptions();
    private static native ByteBuffer getRefContigNames( long indexAddress );
    private static native ByteBuffer createAlignments( ByteBuffer seqs, long indexAddress, ByteBuffer opts );
    static native ByteBuffer createByteBuffer( int size );
    static native void destroyByteBuffer( ByteBuffer alignments );
    private static native String getVersion();
}
