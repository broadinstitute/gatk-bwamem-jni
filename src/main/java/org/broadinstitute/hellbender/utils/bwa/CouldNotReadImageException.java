package org.broadinstitute.hellbender.utils.bwa;

/**
 * Created by valentin on 7/30/17.
 */
public class CouldNotReadImageException extends RuntimeException {

    private final String image;

    public CouldNotReadImageException(final String image, final String message) {
        super(composeMessage(image, message));
        this.image = image;
    }

    public String getImage() {
        return image;
    }

    public CouldNotReadImageException(final String image) {
        super(composeMessage(image, null));
        this.image = image;
    }

    private static String composeMessage(final String image, final String message) {
        return String.format("cannot read the image file '%s'" + (message == null ? "" : ": " + message), image);
    }
}
