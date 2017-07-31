package org.broadinstitute.hellbender.utils.bwa;

/**
 * Created by valentin on 7/30/17.
 */
public class CouldNotReadIndexException extends RuntimeException {

    private final String index;

    public CouldNotReadIndexException(final String index, final String message) {
        super(composeMessage(index, message));
        this.index = index;
    }

    public String getImage() {
        return index;
    }

    public CouldNotReadIndexException(final String index) {
        super(composeMessage(index, null));
        this.index = index;
    }

    private static String composeMessage(final String index, final String message) {
        return String.format("could not read index '%s'" + ((message == null) ? "" : ": " + message), index);
    }
}
