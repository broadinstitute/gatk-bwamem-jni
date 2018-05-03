package org.broadinstitute.hellbender.utils.bwa;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.*;
import java.util.stream.Stream;

@Test
public final class BwaMemIndexTest {
    private static BwaMemIndex index;

    private static final String[] INDEX_EXTENSIONS = {".amb", ".ann", ".bwt", ".pac", ".sa" };

    @BeforeClass
    void openIndex() {
        final String indexImageFile = "src/test/resources/ref.fa.img";
        new File(indexImageFile).deleteOnExit();
        BwaMemIndex.createIndexImageFromIndexFiles("src/test/resources/ref.fa", indexImageFile );
        index = new BwaMemIndex(indexImageFile);
    }

    @AfterClass
    void closeIndex() { index.close(); index = null; }

    @Test
    void testOptsSize() {
        try ( final BwaMemAligner aligner = new BwaMemAligner(index) ) {
            Assert.assertEquals(aligner.getOptsSize(), aligner.getExpectedOptsSize());
        }
    }

    @Test
    void testSimple() {
        final BwaMemAligner aligner = new BwaMemAligner(index);
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(Collections.singletonList(
                // first line in ref.fa
                "GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT".getBytes()
        ));
        Assert.assertNotNull(alignments);
        Assert.assertEquals(alignments.size(), 1);
        final List<BwaMemAlignment> alignmentList = alignments.get(0);
        Assert.assertNotNull(alignmentList);
        Assert.assertEquals(alignmentList.size(), 1);
        testAlignment(alignmentList.get(0), 0, 70, 0, 70, "70M", 0, 0);
    }

    @Test
    void testMulti() {
        final List<String> seqs = new ArrayList<>();
        seqs.add("GGCTTTTAATGCTTTTCAGTGCTAGGTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT"); // 3 snvs
        seqs.add("AATAATAGAGCTTACCATCTGCTGAGTAGACTCCATCTTGAGCAGCAACCACTGAAAAGCATTAAAAGCC"); // rc
        seqs.add("AATACTTCTTTTGAAGCTGCAGTTGTTGCTGCCTTCAACATTAGAATTAATGGGTATTCAATATGATT"); // 2-base deletion
        final BwaMemAligner aligner = new BwaMemAligner(index);
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(seqs,String::getBytes);
        Assert.assertNotNull(alignments);
        Assert.assertEquals(alignments.size(), 3);
        List<BwaMemAlignment> alignmentList = alignments.get(0);
        Assert.assertNotNull(alignmentList);
        Assert.assertEquals(alignmentList.size(), 1);
        testAlignment(alignmentList.get(0), 0, 70, 0, 70, "70M", 3, 0); // 3 snvs
        alignmentList = alignments.get(1);
        Assert.assertNotNull(alignmentList);
        Assert.assertEquals(alignmentList.size(), 1);
        testAlignment(alignmentList.get(0), 0, 70, 0, 70, "70M", 0, 0x10); // rc
        alignmentList = alignments.get(2);
        Assert.assertNotNull(alignmentList);
        Assert.assertEquals(alignmentList.size(), 1);
        testAlignment(alignmentList.get(0), 70, 140, 0, 68, "32M2D36M", 2, 0); // 2-base deletion
    }

    @Test(dataProvider = "testPairData")
    void testPair(final int defaultSetOrClearPEStats) {
        final List<String> seqs = new ArrayList<>();
        seqs.add("GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT"); // ref.fa line 1
        seqs.add("TTGTTTTTAACACCAGAGTCATCCATCACATAATCAAATTTACTTTTAACTCTGGTAAATACTTCATTGT"); // rc ref.fa line 3
        final BwaMemAligner aligner = new BwaMemAligner(index);
        aligner.alignPairs();
        switch (defaultSetOrClearPEStats) {
            case 1:
                aligner.setProperPairEndStats(new BwaMemPairEndStats(200, 10, 1, 600));
                break;
            case 2:
                aligner.dontInferPairEndStats();
                break;
            default:
                aligner.inferPairEndStats();
        }
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(seqs,String::getBytes);
        Assert.assertNotNull(alignments);
        Assert.assertEquals(alignments.size(), 2);
        List<BwaMemAlignment> alignmentList = alignments.get(0);
        Assert.assertNotNull(alignmentList);
        Assert.assertEquals(alignmentList.size(), 1);
        BwaMemAlignment alignment = alignmentList.get(0);
        testAlignment(alignment, 0, 70, 0, 70, "70M", 0, defaultSetOrClearPEStats == 1 ? 0x63 : 0x61);
        Assert.assertEquals(alignment.getMateRefStart(), 140);
        Assert.assertEquals(alignment.getTemplateLen(), 210);
        alignmentList = alignments.get(1);
        Assert.assertNotNull(alignmentList);
        Assert.assertEquals(alignmentList.size(), 1);
        alignment = alignmentList.get(0);
        testAlignment(alignment, 140, 210, 0, 70, "70M", 0, defaultSetOrClearPEStats == 1 ? 0x93 : 0x91);
        Assert.assertEquals(alignment.getMateRefStart(), 0);
        Assert.assertEquals(alignment.getTemplateLen(), -210);
    }

    @DataProvider(name = "testPairData")
    public Object[][] testPairData() {
        final List<Object[]> result = new ArrayList<>(3);
        result.add(new Object[] { 0 });
        result.add(new Object[] { 1 });
        result.add(new Object[] { 2 });
        return result.toArray(new Object[result.size()][]);
    }

    void testAlignment( final BwaMemAlignment alignment,
                        final int refStart, final int refEnd, final int seqStart, final int seqEnd,
                        final String cigar, final int nMismatches, final int samFlag ) {
        Assert.assertEquals(alignment.getRefStart(), refStart);
        Assert.assertEquals(alignment.getRefEnd(), refEnd);
        Assert.assertEquals(alignment.getSeqStart(), seqStart);
        Assert.assertEquals(alignment.getSeqEnd(), seqEnd);
        Assert.assertEquals(alignment.getCigar(), cigar);
        Assert.assertEquals(alignment.getNMismatches(), nMismatches);
        Assert.assertEquals(alignment.getRefId(), 0);
        Assert.assertEquals(alignment.getSamFlag(), samFlag);
    }

    @Test(dataProvider = "algorithmsData")
    void testIndexReference(final BwaMemIndex.Algorithm alg) throws IOException, IOException {
        final char[] refSeq1 = new char[45212];
        final char[] refSeq2 = new char[13415];
        final char[] nucletotides = new char[] { 'A', 'C', 'G', 'T'};
        final Random rdn = new Random(13);
        for (int i = 0; i < refSeq1.length; i++) {
            refSeq1[i] = nucletotides[rdn.nextInt(nucletotides.length)];
        }
        for (int i = 0; i < refSeq2.length; i++) {
            refSeq2[i] = nucletotides[rdn.nextInt(nucletotides.length)];
        }

        final File fastaFile = File.createTempFile("test", ".fasta");
        //fastaFile.deleteOnExit();

        final int basesPerLine = 60;
        try (final PrintWriter fastaWriter = new PrintWriter(new FileWriter(fastaFile))) {
            fastaWriter.println(">seq1");
            for (int i = 0; i < refSeq1.length; i += basesPerLine) {
                fastaWriter.println(new String(refSeq1, i, Math.min(basesPerLine, refSeq1.length - i)));
            }
            fastaWriter.println(">seq2");
            for (int i = 0; i < refSeq2.length; i += basesPerLine) {
                fastaWriter.println(new String(refSeq2, i, Math.min(basesPerLine, refSeq2.length - i)));
            }
        }
        final File imageFile = new File(fastaFile.getPath() + ".idx");
        imageFile.deleteOnExit();
        BwaMemIndex.createIndexImageFromFastaFile(fastaFile.getPath(), imageFile.getPath());
        final BwaMemIndex index = new BwaMemIndex(imageFile.getPath()  );
        Assert.assertEquals(index.getReferenceContigNames(), Arrays.asList("seq1", "seq2"));
        index.close();
        fastaFile.delete();
        imageFile.delete();
    }

    private static void deleteOnExit(final String prefix) {
        Stream.of(INDEX_EXTENSIONS)
                .map(ext -> prefix + ext)
                .map(File::new)
                .forEach(File::deleteOnExit);
    }
    private static void delete(final String prefix) {
        Stream.of(INDEX_EXTENSIONS)
                .map(ext -> prefix + ext)
                .map(File::new)
                .forEach(File::deleteOnExit);
    }

    private static void assertNonEmptyReadableIndexFiles(final String prefix) {
        Stream.of(INDEX_EXTENSIONS)
                .map(ext -> prefix + ext)
                .map(File::new)
                .forEach(f -> {
                    Assert.assertTrue(f.exists());
                    Assert.assertTrue(f.isFile());
                    Assert.assertTrue(f.canRead());
                    Assert.assertTrue(f.length() > 0);
                });
    }

    @DataProvider(name = "algorithmsData")
    Object[][] algorithmsData() {
        return Stream.of(BwaMemIndex.Algorithm.values())
                .map(a -> new Object[] { a }).toArray(Object[][]::new);
    }
}
