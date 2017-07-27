package org.broadinstitute.hellbender.utils.bwa;

/**
 * Created by valentin on 7/30/17.
 */
public class CouldNotCreateIndexImageException extends RuntimeException {

    private final String image;

    public CouldNotCreateIndexImageException(final String image, final String message) {
        super(composeMessage(image, message));
        this.image = image;
    }

    public CouldNotCreateIndexImageException(final String image, final String message, final Throwable cause) {
        super(composeMessage(image, message), cause);
        this.image = image;
    }

    private static String composeMessage(final String image, final String message) {
        return String.format("could not create index image '%s'" + ((message == null) ? "" : ": " + message), image);
    }
}
