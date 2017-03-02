package org.broadinstitute.hellbender.utils.bwa;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

/**
 * Given an open index, this class lets you do alignment of sequences.
 * Don't forget to close it, or you'll leak a little memory.
 * Usage pattern:
 *   Create a BwaMemAligner on some BwaMemIndex
 *   Set your bwa-mem parameters
 *   Align 1 or more chunks of sequences with alignSeqs
 *   Close the BwaMemAligner
 * This class is not thread-safe, but it's very light-weight:  just use a separate instance in each thread.  
 */
public final class BwaMemAligner implements AutoCloseable {
    private final BwaMemIndex index;
    private ByteBuffer opts;

    public BwaMemAligner( final BwaMemIndex index ) {
        this.index = index;
        if ( !index.isOpen() ) {
            throw new IllegalStateException("Can't create aligner: bwa-mem index has been closed");
        }
        opts = BwaMemIndex.createDefaultOptions();
        opts.order(ByteOrder.nativeOrder()).position(0).limit(opts.capacity());
    }

    public boolean isOpen() { return opts != null; }

    @Override
    public void close() {
        if ( opts != null ) {
            BwaMemIndex.destroyByteBuffer(opts);
            opts = null;
        }
    }

    public int getMatchScoreOption() { return getOpts().getInt(0); }
    public void setMatchScoreOption( final int a ) { getOpts().putInt(0, a); }
    public int getMismatchPenaltyOption() { return getOpts().getInt(4); }
    public void setMismatchPenaltyOption( final int b ) { getOpts().putInt(4, b); }
    public int getDGapOpenPenaltyOption() { return getOpts().getInt(8); }
    public void setDGapOpenPenaltyOption( final int o_del ) { getOpts().putInt(8, o_del); }
    public int getDGapExtendPenaltyOption() { return getOpts().getInt(12); }
    public void setDGapExtendPenaltyOption( final int e_del ) { getOpts().putInt(12, e_del); }
    public int getIGapOpenPenaltyOption() { return getOpts().getInt(16); }
    public void setIGapOpenPenaltyOption( final int o_ins ) { getOpts().putInt(16, o_ins); }
    public int getIGapExtendPenaltyOption() { return getOpts().getInt(20); }
    public void setIGapExtendPenaltyOption( final int e_ins ) { getOpts().putInt(20, e_ins); }
    public int getUnpairedPenaltyOption() { return getOpts().getInt(24); }
    public void setUnpairedPenaltyOption( final int pen_unpaired ) { getOpts().putInt(24, pen_unpaired); }
    public int getClip5PenaltyOption() { return getOpts().getInt(28); }
    public void setClip5PenaltyOption( final int pen_clip5 ) { getOpts().putInt(28, pen_clip5); }
    public int getClip3PenaltyOption() { return getOpts().getInt(32); }
    public void setClip3PenaltyOption( final int pen_clip3 ) { getOpts().putInt(32, pen_clip3); }
    public int getBandwidthOption() { return getOpts().getInt(36); }
    public void setBandwidthOption( final int w ) { getOpts().putInt(36, w); }
    public int getZDropOption() { return getOpts().getInt(40); }
    public void setZDropOption( final int zdrop ) { getOpts().putInt(40, zdrop); }
    public long getMaxMemIntvOption() { return getOpts().getLong(48); }
    public void setMaxMemIntvOption( final long max_mem_intv ) { getOpts().putLong(48, max_mem_intv); }
    public int getOutputScoreThresholdOption() { return getOpts().getInt(56); }
    public void setOutputScoreThresholdOption( final int T ) { getOpts().putInt(56, T); }

    public void alignPairs() { setFlagOption(MEM_F_PE|getFlagOption()); }

    // flag bits for the flag option
    public static final int MEM_F_PE = 0x2; // this one's particularly important -- the "paired ends" flag
    public static final int MEM_F_NOPAIRING = 0x4;
    public static final int MEM_F_ALL = 0x8;
    public static final int MEM_F_NO_MULTI = 0x10;
    public static final int MEM_F_NO_RESCUE = 0x20;
    public static final int MEM_F_REF_HDR = 0x100;
    public static final int MEM_F_SOFTCLIP = 0x200;
    public static final int MEM_F_SMARTPE = 0x400;
    public static final int MEM_F_PRIMARY5 = 0x800;
    public int getFlagOption() { return getOpts().getInt(60); }
    public void setFlagOption( final int flag ) { getOpts().putInt(60, flag); }

    public int getMinSeedLengthOption() { return getOpts().getInt(64); }
    public void setMinSeedLengthOption( final int min_seed_len ) { getOpts().putInt(64, min_seed_len); }
    public int getMinChainWeightOption() { return getOpts().getInt(68); }
    public void setMinChainWeightOption( final int min_chain_weight ) { getOpts().putInt(68, min_chain_weight); }
    public int getMaxChainExtendOption() { return getOpts().getInt(72); }
    public void setMaxChainExtendOption( final int max_chain_extend ) { getOpts().putInt(72, max_chain_extend); }
    public float getSplitFactorOption() { return getOpts().getFloat(76); }
    public void setSplitFactorOption( final float split_factor ) { getOpts().putFloat(76, split_factor); }
    public int getSplitWidthOption() { return getOpts().getInt(80); }
    public void setSplitWidthOption( final int split_width ) { getOpts().putInt(80, split_width); }
    public int getMaxSeedOccurencesOption() { return getOpts().getInt(84); }
    public void setMaxSeedOccurencesOption( final int max_occ ) { getOpts().putInt(84, max_occ); }
    public int getMaxChainGapOption() { return getOpts().getInt(88); }
    public void setMaxChainGapOption( final int max_chain_gap ) { getOpts().putInt(88, max_chain_gap); }
    public int getNThreadsOption() { return getOpts().getInt(92); }
    public void setNThreadsOption( final int n_threads ) { getOpts().putInt(92, n_threads); }
    public int getChunkSizeOption() { return getOpts().getInt(96); }
    public void setChunkSizeOption( final int chunk_size ) { getOpts().putInt(96, chunk_size); }
    public float getMaskLevelOption() { return getOpts().getFloat(100); }
    public void setMaxLevelOption( final float max_level ) { getOpts().putFloat(100, max_level); }
    public float getDropRatioOption() { return getOpts().getFloat(104); }
    public void setDropRatioOption( final float drop_ratio ) { getOpts().putFloat(104, drop_ratio); }
    public float getXADropRatio() { return getOpts().getFloat(108); }
    public void setXADropRatio( final float XA_drop_ratio ) { getOpts().putFloat(108, XA_drop_ratio); }
    public float getMaskLevelRedunOption() { return getOpts().getFloat(112); }
    public void setMaskLevelRedunOption( final float max_level_redun ) { getOpts().putFloat(112, max_level_redun); }
    public float getMapQCoefLenOption() { return getOpts().getFloat(116); }
    public void setMapQCoefLenOption( final float mapQ_coef_len ) { getOpts().putFloat(116, mapQ_coef_len); }
    public int getMapQCoefFacOption() { return getOpts().getInt(120); }
    public void setMapQCoefFacOption( final int mapQ_coef_fac ) { getOpts().putInt(120, mapQ_coef_fac); }
    public int getMaxInsOption() { return getOpts().getInt(124); }
    public void setMaxInsOption( final int max_ins ) { getOpts().putInt(124, max_ins); }
    public int getMaxMateSWOption() { return getOpts().getInt(128); }
    public void setMaxMateSWOption( final int max_matesw ) { getOpts().putInt(128, max_matesw); }
    public int getMaxXAHitsOption() { return getOpts().getInt(132); }
    public void setMaxXAHitsOption( final int max_XA_hits ) { getOpts().putInt(132, max_XA_hits); }
    public int getMaxXAHitsAltOption() { return getOpts().getInt(136); }
    public void setMaxXAHitsAltOption( final int max_XA_hits_alt ) { getOpts().putInt(136, max_XA_hits_alt); }
    public byte[] getScoringMatrixOption() {
        final byte[] result = new byte[25];
        final ByteBuffer tmpOpts = getOpts();
        tmpOpts.position(140);
        tmpOpts.get(result);
        return result; }
    public void setScoringMatrixOption( final byte[] mat ) {
        final ByteBuffer tmpOpts = getOpts();
        tmpOpts.position(140);
        tmpOpts.put(mat);
    }
    int getExpectedOptsSize() { return 168; }
    int getOptsSize() { return getOpts().capacity(); }

    public void setIntraCtgOptions() {
        setDGapOpenPenaltyOption(16);
        setIGapOpenPenaltyOption(16);
        setMismatchPenaltyOption(9);
        setClip5PenaltyOption(5);
        setClip3PenaltyOption(5);
    }

    public BwaMemIndex getIndex() {
        return index;
    }

    /**
     * Just align some sequences.
     * @param sequences A list of byte[]'s that contain base calls (ASCII 'A', 'C', 'G', or 'T').
     * @return A list of the same length as the input list.  Each element is a list of alignments for the corresponding sequence.
     */
    public List<List<BwaMemAlignment>> alignSeqs( final List<byte[]> sequences ) {
        return alignSeqs(sequences, seq -> seq);
    }

    /**
     * A more abstract version that takes an iterable of things that can be turned into a byte[] of base calls.
     * @param iterable An iterable over something like a read, that contains a sequence.
     * @param func A lambda that picks the sequence out of your read-like thing.
     * @param <T> The read-like thing.
     * @return A list of (possibly multiple) alignments for each input sequence.
     */
    public <T> List<List<BwaMemAlignment>> alignSeqs( final Iterable<T> iterable, final Function<T,byte[]> func ) {
        final ByteBuffer tmpOpts = getOpts();
        index.refIndex(); // tell the index that we're doing some aligning so that it can't be closed
        final ByteBuffer alignsBuf;
        int nSequences = 0;
        try {
            int bufferCapacity = 4; // buffer will have a 4-byte sequence count as it's first element
            for ( final T ele : iterable ) {
                nSequences += 1;
                bufferCapacity += func.apply(ele).length + 1; // sequence length bytes + 1 for the trailing null
            }
            final ByteBuffer contigBuf = ByteBuffer.allocateDirect(bufferCapacity);
            contigBuf.order(ByteOrder.nativeOrder());
            contigBuf.putInt(nSequences);
            for ( final T ele : iterable ) {
                contigBuf.put(func.apply(ele)).put((byte) 0);
            }
            contigBuf.flip();
            alignsBuf = index.doAlignment(contigBuf, tmpOpts);
        }
        finally {
            index.deRefIndex();
        }
        alignsBuf.order(ByteOrder.nativeOrder()).position(0).limit(alignsBuf.capacity());
        final List<List<BwaMemAlignment>> allAlignments = new ArrayList<>(nSequences);
        while ( nSequences-- > 0 ) {
            int nAligns = alignsBuf.getInt();
            final List<BwaMemAlignment> alignments = new ArrayList<>(nAligns);
            while ( nAligns-- > 0 ) {
                final int flag_mapQ = alignsBuf.getInt();
                final int flags = flag_mapQ >>> 16;
                final int mapQual = flag_mapQ & 0xff;
                final int refId;
                final int refStart;
                final int refEnd;
                final int seqStart;
                final int seqEnd;
                final int nMismatches;
                final int alignerScore;
                final int suboptimalScore;
                final StringBuilder cigar = new StringBuilder();
                final String mdTag;
                final String xaTag;

                // if unmapped
                if ( (flags & 0x4) != 0 ) {
                    refId = -1;
                    refStart = -1;
                    refEnd = -1;
                    seqStart = -1;
                    seqEnd = -1;
                    nMismatches = 0;
                    alignerScore = 0;
                    suboptimalScore = 0;
                    mdTag = null;
                    xaTag = null;
                }
                else { // mapped
                    refId = alignsBuf.getInt();
                    refStart = alignsBuf.getInt();
                    nMismatches = alignsBuf.getInt();
                    alignerScore = alignsBuf.getInt();
                    suboptimalScore = alignsBuf.getInt();
                    int nCigarOps = alignsBuf.getInt();
                    final String cigarOps = "MID?S???????????";
                    if ( nCigarOps <= 0 ) {
                        seqStart = 0;
                        seqEnd = 0;
                        refEnd = refStart;
                    }
                    else {
                        int refLen = 0;
                        int seqLen = 0;
                        int lenOp = alignsBuf.getInt();
                        int len = lenOp >>> 4;
                        char op = cigarOps.charAt(lenOp & 0x0f);
                        cigar.append(len);
                        cigar.append(op);
                        seqStart = op == 'S' ? len : 0;
                        if ( op == 'M' || op == 'D' ) refLen += len;
                        if ( op == 'M' || op == 'I' ) seqLen += len;
                        while ( --nCigarOps > 0 ) {
                            lenOp = alignsBuf.getInt();
                            len = lenOp >>> 4;
                            op = cigarOps.charAt(lenOp & 0x0f);
                            cigar.append(len);
                            cigar.append(op);
                            if ( op == 'M' || op == 'D' ) refLen += len;
                            if ( op == 'M' || op == 'I' ) seqLen += len;
                        }
                        refEnd = refStart + refLen;
                        seqEnd = seqStart + seqLen;
                    }
                    mdTag = getTag(alignsBuf);
                    xaTag = getTag(alignsBuf);
                }

                final int mateRefId;
                final int mateStartPos;
                final int templateLen;

                // if unpaired, or mate unmapped
                if ( (flags & 0x1) == 0 || (flags & 0x8) != 0 ) {
                    mateRefId = -1;
                    mateStartPos = -1;
                    templateLen = 0;
                } else { // has mapped mate
                    mateRefId = alignsBuf.getInt();
                    mateStartPos = alignsBuf.getInt();
                    templateLen = alignsBuf.getInt();
                }
                alignments.add(new BwaMemAlignment(flags, refId, refStart, refEnd, seqStart, seqEnd, mapQual,
                        nMismatches, alignerScore, suboptimalScore, cigar.toString(), mdTag, xaTag,
                        mateRefId, mateStartPos, templateLen));
            }
            allAlignments.add(alignments);
        }
        BwaMemIndex.destroyByteBuffer(alignsBuf);
        return allAlignments;
    }

    private String getTag( final ByteBuffer buffer ) {
        int tagLen = buffer.getInt();
        if ( tagLen == 0 ) return null;
        byte[] tagBytes = new byte[(tagLen+3)&~3];
        buffer.get(tagBytes);
        return new String(tagBytes, 0, tagLen);
    }

    private ByteBuffer getOpts() {
        if ( opts == null ) {
            throw new IllegalStateException("The aligner has been closed.");
        }
        return opts;
    }
}
