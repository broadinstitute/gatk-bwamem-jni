package org.broadinstitute.hellbender.utils.bwa;

public class BwaMemException extends RuntimeException {
    public BwaMemException( final String message ) { super(message); }
    public BwaMemException( final String message, final Exception cause ) { super(message, cause); }
}
