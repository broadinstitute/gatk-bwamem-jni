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
public final class BwaMemPairEndStats {

    /**
     * Default number of std. deviation from the mean for the lowest and highest insert size that is considered normal.
     * <p>
     *     Its value, {@value #DEFAULT_LOW_AND_HIGH_SIGMA}, corresponds to the one use in bwa-mem source code at the
     *     time this code was written.
     * </p>
     */
    public static final int DEFAULT_LOW_AND_HIGH_SIGMA = 4;

    /**
     * Default std. deviation to average ratio to use in case the deviation is not provided.
     * <p>
     *     Its value, {@value #DEFAULT_STD_TO_AVERAGE_RATIO}, corresponds to the one use in bwa-mem source code at the
     *     time this code was written.
     * </p>
     */
    public static final double DEFAULT_STD_TO_AVERAGE_RATIO = .1;

    /**
     * Constant to indicate that pair-end insert size inference failed or that it should not be inferred depending on
     * the context.
     */
    public static final BwaMemPairEndStats FAILED = new BwaMemPairEndStats();

    /**
     * Constant to indicate that pair-end insert size inference failed or hat is should not be inferred depending on
     * the context.
     */
    public static final BwaMemPairEndStats DO_NOT_INFER = FAILED;

    /**
     * The shortest insert size that is considered normal.
     * <br/>
     * Only defined if {@link #failed} is {@code false}, and if so is always in the [1..{@link #average}] range.
     */
    public final int low;

    /**
     * The longest insert size that is considered normal.
     * <br/>
     * Only defined if {@link #failed} is {@code false}, and if so is always in the [#average .. +Inf) range.
     */
    public final int high;

    /**
     * Depending on the context, indicates iff the insert-size inference failed or
     * if we are requesting not to perform such a inference.
     * <br/>
     * When this field is {@code true}, nothing can be expected on what values other fields take.
     */
    public final boolean failed;

    /**
     * The average insert size.
     * <br/>
     * Only defined if {@link #failed} is {@code false}, and if so is always in the range [1 .. +Inf).
     */
    public final double average;

    /**
     * The insert size std. deviation.
     * <br/>
     * Only defined if {@link #failed} is {@code false}, and if so is finite and equal or greater than 0.
     */
    public final double std;


    /**
     * Composes a new instance given the estimated average and std. deviation
     * <p>
     *     low and high limits are calculated using the recipe use in bwa mem's code-base
     *     where the 4th sigma is used as the upper limit rounding it to the closest integer.
     * </p>
     * <p>
     *     The resulting stats object's {@link #failed} field will be set to {@code false}.
     * </p>
     * @param average the insert size average estimate.
     * @param std the insert size standard deviation estimate.
     * @throws IllegalArgumentException if either input is {@link Double#NaN}, infinite. Also if the average is
     *  less than 1.0 and if the std. deviation is strictly negative.
     */
    public BwaMemPairEndStats(final double average, final double std) {
        this(average, std,
         /* low  = */ Math.max(1, (int) Math.round(average - DEFAULT_LOW_AND_HIGH_SIGMA * std)),
         /* high = */ Math.max(1, (int) Math.round(average + DEFAULT_LOW_AND_HIGH_SIGMA * std)));
    }

    /**
     * Composes a new instance given the estimated average.
     * <p>
     * The std. deviation would be calculated out of this mean using the default deviation
     * to average ratio {@link #DEFAULT_STD_TO_AVERAGE_RATIO}.
     * </p>
     * <p>
     *     The resulting stats object's {@link #failed} field will be set to {@code false}.
     * </p>
     * <p>
     *     Low and high insert-sizes will be estimated using the approach described in {@link #BwaMemPairEndStats(double, double)}.
     * </p>
     * <p>
     *     The resulting stats object's {@link #failed} field will be set to {@code false}.
     * </p>
     * @param average the insert size average estimate.
     * @throws IllegalArgumentException if the input average is {@link Double#NaN}, infinite or less than 1.
     */
    public BwaMemPairEndStats(final double average) {
        this(average, average * DEFAULT_STD_TO_AVERAGE_RATIO);
    }

    /**
     * Composes a new instance giving arbitrary (albeit valid and consistent) values to all its numerical
     * fields: {@link #average}, {@link #std}, {@link #low} and {@link #high}.
     * <p>
     *     The resulting stats object's {@link #failed} field will be set to {@code false}.
     * </p>
     * @param average the insert-size average in [1 .. +Inf) range.
     * @param std the insert-size std. dev in [0 .. +Inf) range.
     * @param low the shortest insert-size to be considered "normal" [1 .. {@code average}]
     * @param high the longest insert-size to be considered "normal" [{@code average} .. +Inf).
     *
     * @throws IllegalArgumentException if the values provided are invalid or inconsistent based on the constraints
     * listed above.
     */
    public BwaMemPairEndStats(final double average, final double std, final int low, final int high) {
        if (Double.isNaN(average) || !Double.isFinite(average) || average < 1) {
            throw new IllegalArgumentException("invalid input average: " + average);
        } else if (Double.isNaN(std) || !Double.isFinite(std) || std < 0) {
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

    /**
     * Constructor used to create the singleton {@link #FAILED} instance.
     */
    private BwaMemPairEndStats() {
        this.failed = true;
        this.average = Double.NaN;
        this.std = Double.NaN;
        this.low = Integer.MAX_VALUE;
        this.high = Integer.MIN_VALUE;
    }

    @Override
    public String toString() {
        if (failed) {
            return "InsertSize ~ FAILED/DO_NOT_INFER";
        }
        return String.format("InsertSize ~ N(%.2f, %.2f) in [%d, %d]", average, std, low, high);
    }

    @Override
    public boolean equals(final Object other) {
        return other instanceof BwaMemPairEndStats && equals((BwaMemPairEndStats) other);
    }

    @Override
    public int hashCode() {
        if (this.failed) {
            return 0;
        } else {
            return ((((Double.hashCode(average) * 31) +
                    Double.hashCode(std) * 31) +
                    Integer.hashCode(low) * 31) +
                    Integer.hashCode(high) * 31);
        }
    }

    public boolean equals(final BwaMemPairEndStats other) {
        if (other == null || this.failed != other.failed) {
            return false;
        } else {
            return this.failed || (this.average == other.average &&
                    this.std == other.std &&
                    this.low == other.low && this.high == other.high);
        }
    }


}
