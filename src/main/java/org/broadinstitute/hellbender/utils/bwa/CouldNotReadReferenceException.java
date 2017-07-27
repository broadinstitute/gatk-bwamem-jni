package org.broadinstitute.hellbender.utils.bwa;

/**
 * Created by valentin on 7/30/17.
 */
public class CouldNotReadReferenceException extends RuntimeException {

    private final String reference;

    public CouldNotReadReferenceException(final String reference, final String message) {
        super(composeMessage(reference, message));
        this.reference = reference;
    }

    public CouldNotReadReferenceException(final String reference) {
        super(composeMessage(reference, null));
        this.reference = reference;
    }

    public String getReference() {
        return reference;
    }

    private static String composeMessage(final String reference, final String message) {
        return String.format("cannot read the reference file '%s'" + (message == null ? "" : ": " + message), reference);
    }
}
