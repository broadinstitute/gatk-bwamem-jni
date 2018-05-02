package org.broadinstitute.hellbender.utils.bwa;

/**
 * Pair-end read orientations where "F" stands for forward and "R" for reverse.
 */
public enum BwaMemPairEndOrientation {

    /**
     * Forward-Foward.
     * <br/>
     * Both mates align on the reference strand. Possible inversion when dealing with
     * Illumina short-read data.
     */
    FF,

    /**
     * Forward-Reverse.
     * <br/>
     * First mate (upstream) maps to the foward strand and its mate (downstream). This
     * the common case in Illumina short-read data in the absence of structural variation.
     */
    FR,

    /**
     * Reverse-Forward.
     * <br/>
     * Both mates look away from each other; the upstream pair in mapped to the reverse, and the downstream
     * mate is mapped to the forward strand. In Illumina data this is considered and indication of a
     * possible tandem-duplication.
     */
    RF,

    /**
     * Reverse-Reverse
     * <br/>
     * Similar to {@link ##FF} this might indicate the presence of an inversion.
     */
    RR
}
