package org.broadinstitute.hellbender.utils.bwa;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

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

    private static final int MAXIMUM_NUMBER_OF_CHARACTER_BEFORE_FIRST_HEADER_IN_FASTA_FILES =  4092;
    private static final char FASTA_HEADER_PREFIX_CHAR = '>';

    public static final List<String> INDEX_FILE_EXTENSIONS =
            Collections.unmodifiableList(Arrays.asList(
                    ".amb", ".ann", ".bwt", ".pac", ".sa" ));

    public static final String IMAGE_FILE_EXTENSION = ".img";

    public static final List<String> FASTA_FILE_EXTENSIONS =
            Collections.unmodifiableList(Arrays.asList(".fasta", ".fa"));

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
          * Linear-time algorithm for constructing suffix array.
          * It requires 5.37N memory where N is the size of the database.
          * IS is moderately fast, but does not work with database larger than 2GB.
          * IS is the default algorithm due to its simplicity. The current codes for
          * IS algorithm are reimplemented by Yuta Mori.
          *
          * @see "http://bio-bwa.sourceforge.net/bwa.shtml"
          */
        IS,

        /**
         * The <i>Ropebwt2</i> algorithm.
         *
         * @see "https://arxiv.org/pdf/1406.0426.pdf"
         */
        RB2;

        /**
         * The string name use by the Bwa command line to denote that algorithm.
         * @return never {@code null} and unique across algorithms.
         */
        public String toBwaName() {
            return this.toString().toLowerCase();
        }
    }

    private final String indexImageFile; // stash this for error messages
    private volatile long indexAddress; // address where the index was memory-mapped (for use by C code)
    private final AtomicInteger refCount; // keep track of how many threads are actively aligning
    private final List<String> refContigNames; // the reference dictionary from the index
    private static volatile boolean nativeLibLoaded = false; // whether we've loaded the native library or not

    private static String resolveFastaFileExtension(final String fasta) {
        final Optional<String> extension = FASTA_FILE_EXTENSIONS.stream()
                .filter(fasta::endsWith).findFirst();
        if (!extension.isPresent()) {
            throw new IllegalArgumentException(
                    String.format("the fasta file provided '%s' does not have any of the standard fasta extensions: %s",
                            fasta, FASTA_FILE_EXTENSIONS.stream().collect(Collectors.joining(", "))));
        } else {
            return extension.get();
        }
    }

    /**
     * Create the index image file for a complete set of BWA index files
     * @param indexPrefix the location of the index files.
     * @param imageFile the location of the new index image file.
     *
     * <p>
     *     <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt process
     *     interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     *
     * @throws IllegalArgumentException if {@code indexPrefix} is {@code null}
     *  or it does not look like it points to a complete set of index files.
     * @throws IllegalArgumentException if {@code imageFile} is {@code null}.
     */
    public static void createIndexImageFromIndexFiles(final String indexPrefix, final String imageFile) {
        if (indexPrefix == null) {
            throw new IllegalArgumentException("the index prefix cannot be null");
        } else if (imageFile == null) {
            throw new IllegalArgumentException("the image file cannot be null");
        }
        assertLooksLikeIndexPrefix(indexPrefix);
        loadNativeLibrary();
        createIndexImageFile(indexPrefix, imageFile);
    }

    /**
     * Create the index image file for a fasta file.
     * <p>The name of the output index file will be determined by the name of the fasta file using
     * {@link #getDefaultIndexImageNameFromFastaFile(String)}.
     *
     * <p>
     *     This is equivalent to calling
     *     {@code {@link #createIndexImageFromFastaFile(String,String) createIndexImageFromFastaFile(X, getDefaultIndexImageNameFromFastaFile(X)))}}
     * </p>
     * <p>
     *     <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt process
     *     interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     *
     * @param fasta the location of the fasta reference file.
     * @return location of the generated index image file (equivaled to {@link #getDefaultIndexImageNameFromFastaFile(String)}).
     * @throws IllegalArgumentException if {@code fasta} is {@code null} or
     *  does not finish with out of the standard fasta extension names (listed in {@link #FASTA_FILE_EXTENSIONS}).
     * @throws InvalidFileFormatException if {@code fasta} does not seem to be a fasta formatted regular readable file.
     * @throws CouldNotCreateIndexImageException if for some reason we could not create the index file.
     * @see #getDefaultIndexImageNameFromFastaFile(String)
     */
    public static String createIndexImageFromFastaFile(final String fasta) {
        final String imageFile = getDefaultIndexImageNameFromFastaFile(fasta);
        createIndexImageFromFastaFile(fasta, imageFile);
        return imageFile;
    }

    /**
     * Gets the default index image name for the provided FASTA.
     * <p>The default name is substitutes its extension (typically <i>.fasta</i> or <i>.fa</i>)
     * by {@link #IMAGE_FILE_EXTENSION} (<i>.img</i>).</p>
     *
     *
     * @param fasta the location of the fasta reference file.
     * @throws IllegalArgumentException if {@code fasta} is {@code null} or
     *  does not finish with out of the standard fasta extension names (listed in {@link #FASTA_FILE_EXTENSIONS}).
     */
    public static String getDefaultIndexImageNameFromFastaFile(final String fasta) {
        if (fasta == null) {
            throw new IllegalArgumentException("the input fasta file name cannot be null");
        }
        final String extension = resolveFastaFileExtension(fasta);
        final String prefix = fasta.substring(0, fasta.length() - extension.length());
        return prefix + IMAGE_FILE_EXTENSION;
    }

    /**
     * Creates the index image file for a reference fasta file.
     * <p>
     *     The index will be created using the default algorithm {@link Algorithm#AUTO}.
     * </p>
     * <p>
     *     Calling this method is equivalent to calling {@link #createIndexImageFromFastaFile(String, String, Algorithm) createIndexImageFromFastaFile(a, b, Algorithm.AUTO)}
     * <p>
     *     <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt process
     *     interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     * @param fasta the location of the targeted reference.
     * @param imageFile the location of the new index image file.
     * @throws IllegalArgumentException if {@code fasta} is {@code null}
     *  or it does not look like it points to a fasta formatted readable file.
     * @throws IllegalArgumentException if {@code imageFile} is {@code null}.
     * @throws InvalidFileFormatException if {@code fasta} does not seem to be
     * a fasta formatted regular and readable file.
     * @throws CouldNotCreateIndexImageException if there was a problem creating
     * the output image.
     * @throws CouldNotCreateIndexException if there was some problem while creating
     * the intermediary index file set.
     */
    public static void createIndexImageFromFastaFile(final String fasta, final String imageFile) {
        createIndexImageFromFastaFile(fasta, imageFile, Algorithm.AUTO);
    }

    /**
     * Creates the index image file for a reference fasta file.
     * <p>
     *     <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt process
     *     interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     * @param fasta the location of the targeted reference.
     * @param imageFile the location of the new index image file.
     * @param algo the algorithm to use to construct the index (see {@link Algorithm} to see what there is available.).
     * @throws IllegalArgumentException if {@code fasta} is {@code null}
     *  or it does not look like it points to a fasta formatted readable file.
     * @throws IllegalArgumentException if {@code imageFile} is {@code null}.
     * @throws InvalidFileFormatException if {@code fasta} does not seem to be
     * a fasta formatted regular and readable file.
     * @throws CouldNotCreateIndexImageException if there was a problem creating
     * the output image.
     * @throws CouldNotCreateIndexException if there was some problem while creating
     * the intermediary index file set.
     */
    public static void createIndexImageFromFastaFile( final String fasta, final String imageFile, final Algorithm algo) {
        assertLooksLikeFastaFile(fasta);
        assertCanCreateOrOverwriteImageFile(imageFile);
        if (algo == null) {
            throw new IllegalArgumentException("the input algorithm must not be null");
        }

        final File indexPrefix = createTempIndexPrefix(fasta);
        loadNativeLibrary();
        createReferenceIndex(fasta, indexPrefix.getPath(), algo.toBwaName());
        createIndexImageFile(indexPrefix.getPath(), imageFile);
        deleteIndexFiles(indexPrefix);
    }

    private static void assertCanCreateOrOverwriteImageFile(final String imageFile) {
        if (imageFile == null) {
            throw new IllegalArgumentException("the image file cannot be null");
        } else {
            final File file = new File(imageFile);
            try {
                if (!file.createNewFile()) {
                    if (!file.isFile() || !file.canWrite()) {
                        throw new CouldNotCreateIndexImageException(imageFile, "already exists as a non-regular or unwritable file");
                    }
                } else {
                    file.delete();
                }
            } catch (final IOException ex) {
                throw new CouldNotCreateIndexImageException(imageFile, ex.getMessage(), ex);
            }
        }
    }

    /**
     * Checks whether the input index prefix seems to point to a complete set
     * of readable index files.
     *
     * @param indexPrefix the target index prefix.
     * @throws IllegalArgumentException if {@code indexPrefix} is {@code null}.
     * @throws CouldNotReadIndexException if that is not the case.
     */
    private static void assertLooksLikeIndexPrefix(final String indexPrefix) {
        if (indexPrefix == null) {
            throw new IllegalArgumentException("the input index prefix cannot be null");
        }
        INDEX_FILE_EXTENSIONS.stream()
                .map(ext -> indexPrefix + ext)
                .forEach(file -> assertNonEmptyReadableIndexFile(indexPrefix, file));
    }

    private static void deleteIndexFiles(final File indexPrefix) {
        INDEX_FILE_EXTENSIONS.stream()
                .map(ext -> new File(indexPrefix + ext))
                .forEach(File::delete);
        indexPrefix.delete();
    }

    private static File createTempIndexPrefix(final String fasta) {
        final File indexPrefix;
        try {
            indexPrefix = File.createTempFile("temporal-index","");
        } catch (final IOException ex) {
            throw new CouldNotCreateIndexException(fasta, "no-location","failure to create a temporal file");
        }
        indexPrefix.deleteOnExit();
        INDEX_FILE_EXTENSIONS.stream()
                .map(ext -> new File(indexPrefix + ext))
                .forEach(File::deleteOnExit);
        return indexPrefix;
    }

    private static void assertLooksLikeFastaFile(final String fasta) {
        resolveFastaFileExtension(fasta);
        if (!nonEmptyReadableFile(fasta)) {
            throw new CouldNotReadReferenceException(fasta, "input file unreachable or not a file");
        }
        try (final BufferedReader reader = new BufferedReader(new FileReader(fasta))) {
            int c;
            int offset = 0;
            while (offset++ < MAXIMUM_NUMBER_OF_CHARACTER_BEFORE_FIRST_HEADER_IN_FASTA_FILES && (c = reader.read()) != -1) {
                if (!Character.isSpaceChar(c)) {
                    if (c == FASTA_HEADER_PREFIX_CHAR) {
                        break;
                    } else {
                        throw new InvalidFileFormatException(fasta, "the file provided does not seem to be a fasta file (first non-space character in the first 4K is not '" + FASTA_HEADER_PREFIX_CHAR + "'");
                    }
                }
            }

        } catch (final IOException ex) {
            throw new InvalidFileFormatException(fasta, "problems reading the content of the reference fasta file'", ex);
        }
    }

    /**
     * Loads an index from an image file.
     * <p>
     *     You can use other methods to create such
     *     indexes from fasta reference ({@link #createIndexImageFromFastaFile} or their index files
     *     ({@link #createIndexImageFromIndexFiles}).
     * </p>
     * <p>
     *     <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt process
     *     interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     *
     * @throws IllegalArgumentException if {@code indexImageFile} is {@code null}.
     * @throws CouldNotReadImageException if some problem occurred when loading the
     *  image file.
     */
    public BwaMemIndex( final String indexImageFile ) {
        this.indexImageFile = indexImageFile;
        loadNativeLibrary();
        assertNonEmptyReadableImageFile(indexImageFile);
        refCount = new AtomicInteger();
        indexAddress = openIndex(indexImageFile);
        if ( indexAddress == 0L ) {
            throw new CouldNotReadImageException(indexImageFile, "unable to open bwa-mem index");
        }
        ByteBuffer refContigNamesBuf = getRefContigNames(indexAddress);
        if ( refContigNamesBuf == null ) {
            throw new CouldNotReadImageException("unable to retrieve reference contig names from bwa-mem index");
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

    private void assertNonEmptyReadableImageFile(final String image) {
        if (!nonEmptyReadableFile(image)) {
            throw new CouldNotReadImageException(image, "is empty or is not readable");
        }
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

    /**
     * Close the index and release the (non-Java) memory that's been allocated
     *
     * <p>
     *     <b><i>WARNING!</i></b>: Notice that currently this method is making JNI call that might result in an abrupt process
     *     interruption (e.g. exit or abort system call) and so the control may never be returned.
     * </p>
     */
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
                    destroyIndex(addr);
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

    ByteBuffer doAlignment( final ByteBuffer seqs, final ByteBuffer opts, final BwaMemPairEndStats peStats) {
        final ByteBuffer alignments = createAlignments(seqs, indexAddress, opts, peStats);
        if ( alignments == null ) {
            throw new IllegalStateException("Unable to get alignments from bwa-mem index "+indexImageFile+": We don't know why.");
        }
        return alignments;
    }

    private static void assertNonEmptyReadableIndexFile(final String index, final String fileName ) {
        if ( !nonEmptyReadableFile(fileName) )
            throw new CouldNotReadIndexException(index, "Missing bwa index file: "+ fileName);
    }

    private static boolean nonEmptyReadableFile( final String fileName ) {
        if (fileName == null) {
            throw new IllegalArgumentException("the input file name cannot be null");
        }
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

    private static native boolean createReferenceIndex(String referenceName, String indexPrefix, String algorithmName);
    private static native boolean createIndexImageFile(String indexPrefix, String imageName );
    private static native long openIndex( String indexImageFile );
    private static native int destroyIndex( long indexAddress );
    static native ByteBuffer createDefaultOptions();
    private static native ByteBuffer getRefContigNames( long indexAddress );
    private static native ByteBuffer createAlignments( ByteBuffer seqs, long indexAddress, ByteBuffer opts, BwaMemPairEndStats peStats);
    static native ByteBuffer createByteBuffer( int size );
    static native void destroyByteBuffer( ByteBuffer alignments );
    private static native String getVersion();
}
