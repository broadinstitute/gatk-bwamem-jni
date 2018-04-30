package org.broadinstitute.hellbender.utils.bwa;

/**
 * Equivalent to Bwa's mem_pestat_t struct.
 * At the time this was written such a type declaration looked like this:
 * <pre>
 *    // bwamem.h
 *    typedef struct {
            int low, high;   // lower and upper bounds within which a read pair is considered to be properly paired
            int failed;      // non-zero if the orientation is not supported by sufficient data
            double avg, std; // mean and stddev of the insert size distribution
     } mem_pestat_t;
 * </pre>
 */
public class BwaMemPairEndStats {

    public static final int DEFAULT_LOW_AND_HIGH_SIGMA = 4;
    public static final double DEFAULT_STD_TO_AVERAGE_RATIO = .1;
    public static final BwaMemPairEndStats FAILED = new BwaMemPairEndStats();

    public final int low;
    public final int high;
    public final boolean failed;
    public final double average;
    public final double std;


    /**
     * Composes a new instance given the estimated average and std. deviation
     * <p>
     *     low and high limits are calculated using the recipe use in bwa mem's code-base
     *     where the 4th sigma is used as the upper limit rounding it to the closest integer.
     * </p>
     *
     * @param average
     * @param std
     */
    public BwaMemPairEndStats(final double average, final double std) {
        if (!Double.isFinite(average) || average <= 0) {
            throw new IllegalArgumentException("invalid input average: " + average);
        } else if (!Double.isFinite(std) && std != 0) {
            throw new IllegalArgumentException("invalid std. err: " + std);
        } else {
            final double absStd = std >= 0 ? std : -std;
            this.failed = false;
            this.average = average;
            this.std = absStd;
            this.low = Math.max(1, (int) Math.round(average - DEFAULT_LOW_AND_HIGH_SIGMA * absStd));
            this.high = Math.max(1, (int) Math.round(average + DEFAULT_LOW_AND_HIGH_SIGMA * absStd));
        }
    }

    public BwaMemPairEndStats(final double average) {
        if (!Double.isFinite(average) || average <= 0) {
            throw new IllegalArgumentException("invalid input average: " + average);
        } else {
            this.average = average;
            this.failed = false;
            this.std = this.average * DEFAULT_STD_TO_AVERAGE_RATIO;
            this.low = Math.max(1, (int) Math.round(average - DEFAULT_LOW_AND_HIGH_SIGMA * this.std));
            this.high = Math.max(1, (int) Math.round(average + DEFAULT_LOW_AND_HIGH_SIGMA * this.std));
        }
    }

    public BwaMemPairEndStats(final double average, final double std, final int low, final int high) {
        if (!Double.isFinite(average) || average <= 0) {
            throw new IllegalArgumentException("invalid input average: " + average);
        } else if (!Double.isFinite(std) && std != 0) {
            throw new IllegalArgumentException("invalid std. err: " + std);
        } else if (low > average) {
            throw new IllegalArgumentException("the low limit cannot be larger than the average");
        } else if (high < average) {
            throw new IllegalArgumentException("the high limit cannot be larger than the average");
        } else {
            this.failed = false;
            this.average = average;
            this.std = std;
            this.high = high;
            this.low = low;
        }
    }

    private BwaMemPairEndStats() {
        this.failed = true;
        this.average = Double.NaN;
        this.std = Double.NaN;
        this.low = Integer.MAX_VALUE;
        this.high = Integer.MIN_VALUE;
    }
}
