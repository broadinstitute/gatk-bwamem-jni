package org.broadinstitute.hellbender.utils.bwa;

/**
 * Info from the Aligner about an alignment to reference that it discovered for some sequence.
 * Please note that the refId is with respect to the BWA index reference names.  These won't necessarily agree with
 * reference IDs in SAM or BAM headers.
 */
public class BwaMemAlignment {
    private final int samFlag;     // flag bits per SAM format standard
    private final int refId;       // index into reference dictionary (-1 if unmapped)
    private final int refStart;    // 0-based coordinate, inclusive (-1 if unmapped)
    private final int refEnd;      // 0-based coordinate, exclusive (-1 if unmapped)
    private final int seqStart;    // 0-based coordinate, inclusive (-1 if unmapped)
    private final int seqEnd;      // 0-based coordinate, exclusive (-1 if unmapped)
    private final int mapQual;     // phred-scaled mapping quality (-1 if unmapped)
    private final int nMismatches; // number of mismatches (i.e., value of the NM tag in a SAM/BAM) (-1 if unmapped)
    private final int alignerScore; // for AS tag
    private final int suboptimalScore; // for bwa-specific XS tag
    private final String cigar;     // cigar for alignment (empty if unmapped)
    private final String mdTag;    // the MD tag
    private final String xaTag;    // the XA tag
    private final int mateRefId;   // mate's refId (-1 if unpaired or if mate unmapped)
    private final int mateRefStart;// mate's reference start (-1 if unpaired or if mate unmapped)
    private final int templateLen; // inferred template length (0 if unpaired, mate unmapped, or on different ref contigs)

    public BwaMemAlignment(final int samFlag, final int refId, final int refStart, final int refEnd,
                           final int seqStart, final int seqEnd, final int mapQual,
                           final int nMismatches, final int alignerScore, final int suboptimalScore,
                           final String cigar, final String mdTag, final String xaTag,
                           final int mateRefId, final int mateRefStart, final int templateLen ) {
        this.samFlag = samFlag;
        this.refId = refId;
        this.refStart = refStart;
        this.refEnd = refEnd;
        this.seqStart = seqStart;
        this.seqEnd = seqEnd;
        this.mapQual = mapQual;
        this.nMismatches = nMismatches;
        this.alignerScore = alignerScore;
        this.suboptimalScore = suboptimalScore;
        this.cigar = cigar;
        this.mdTag = mdTag;
        this.xaTag = xaTag;
        this.mateRefId = mateRefId;
        this.mateRefStart = mateRefStart;
        this.templateLen = templateLen;
    }

    public int getSamFlag() { return samFlag; }
    public int getRefId() { return refId; }
    public int getRefStart() { return refStart; }
    public int getRefEnd() { return refEnd; }
    public int getSeqStart() { return seqStart; }
    public int getSeqEnd() { return seqEnd; }
    public int getMapQual() { return mapQual; }
    public int getNMismatches() { return nMismatches; }
    public int getAlignerScore() { return alignerScore; }
    public int getSuboptimalScore() { return suboptimalScore; }
    public String getCigar() { return cigar; }
    public String getMDTag() { return mdTag; }
    public String getXATag() { return xaTag; }
    public int getMateRefId() { return mateRefId; }
    public int getMateRefStart() { return mateRefStart; }
    public int getTemplateLen() { return templateLen; }
}
