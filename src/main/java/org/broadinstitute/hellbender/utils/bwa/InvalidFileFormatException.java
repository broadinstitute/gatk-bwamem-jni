package org.broadinstitute.hellbender.utils.bwa;

import java.util.OptionalLong;

/**
 * Created by valentin on 7/28/17.
 */
public class InvalidFileFormatException extends InvalidInputException {

    private final String file;
    private final OptionalLong line;

    public InvalidFileFormatException(final String file, final OptionalLong line, final String message, final Throwable cause) {
        super(composeMessage(file, line, message), cause);
        this.file = file;
        this.line = line;
    }

    public InvalidFileFormatException(final String file, final String message, final Throwable cause) {
        this(file, OptionalLong.empty(), message, cause);
    }

    public InvalidFileFormatException(final String file, final String message) {
        this(file, OptionalLong.empty(), null, null);
    }

    private static String composeMessage(final String file, final OptionalLong line, final String message) {
        if (file == null) {
            throw new IllegalArgumentException("the input file cannot be null");
        } else if (line.isPresent() && line.getAsLong() < 1) {
            throw new IllegalArgumentException("the line number cannot be null");
        }
        final String location = "file " + file + (line.isPresent() ? (" (" + line.getAsLong() + ")") : "");
        final String details = message == null ? ": invalid format" : ": " + message;
        return location + details;
    }

    public OptionalLong getLine() {
        return line;
    }

    public String getFile() {
        return file;
    }
}
