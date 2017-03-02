package org.broadinstitute.hellbender.utils.bwa;

import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

@Test
public final class BwaMemIndexTest {
    private static BwaMemIndex index;

    @BeforeClass
    void openIndex() {
        final String indexImageFile = "src/test/resources/ref.fa.img";
        new File(indexImageFile).deleteOnExit();
        BwaMemIndex.createIndexImage("src/test/resources/ref.fa", indexImageFile );
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

    @Test
    void testPair() {
        final List<String> seqs = new ArrayList<>();
        seqs.add("GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT"); // ref.fa line 1
        seqs.add("TTGTTTTTAACACCAGAGTCATCCATCACATAATCAAATTTACTTTTAACTCTGGTAAATACTTCATTGT"); // rc ref.fa line 3
        final BwaMemAligner aligner = new BwaMemAligner(index);
        aligner.alignPairs();
        final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(seqs,String::getBytes);
        Assert.assertNotNull(alignments);
        Assert.assertEquals(alignments.size(), 2);
        List<BwaMemAlignment> alignmentList = alignments.get(0);
        Assert.assertNotNull(alignmentList);
        Assert.assertEquals(alignmentList.size(), 1);
        BwaMemAlignment alignment = alignmentList.get(0);
        testAlignment(alignment, 0, 70, 0, 70, "70M", 0, 0x61);
        Assert.assertEquals(alignment.getMateRefStart(), 140);
        Assert.assertEquals(alignment.getTemplateLen(), 210);
        alignmentList = alignments.get(1);
        Assert.assertNotNull(alignmentList);
        Assert.assertEquals(alignmentList.size(), 1);
        alignment = alignmentList.get(0);
        testAlignment(alignment, 140, 210, 0, 70, "70M", 0, 0x91);
        Assert.assertEquals(alignment.getMateRefStart(), 0);
        Assert.assertEquals(alignment.getTemplateLen(), -210);
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
}
