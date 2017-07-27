package org.broadinstitute.hellbender.utils.bwa;

/**
 * Created by valentin on 7/29/17.
 */
public class CouldNotCreateIndexException extends RuntimeException {
    public CouldNotCreateIndexException(final String fasta, final String index, final String message) {
        super(composeMessage(fasta, index, message));
    }

    private static String composeMessage(final String fasta, final String prefix, final String message) {
        return String.format("Failed to create index for index file '%s' in location '%s'" + message == null ? "" : ": " + message, fasta, prefix);
    }
}
