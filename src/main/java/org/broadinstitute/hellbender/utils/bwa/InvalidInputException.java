package org.broadinstitute.hellbender.utils.bwa;

/**
 * Created by valentin on 7/28/17.
 */
public class InvalidInputException extends RuntimeException {

    public InvalidInputException() {
    }

    public InvalidInputException(final String message) {
        super(message);
    }

    public InvalidInputException(final String message, final Throwable cause) {
        super(message, cause);
    }

    public InvalidInputException(final Throwable cause) {
        super(cause);
    }

    public InvalidInputException(final String message, final Throwable cause, final boolean enableSuppression, boolean writableStackTrace) {
        super(message, cause, enableSuppression, writableStackTrace);
    }
}
