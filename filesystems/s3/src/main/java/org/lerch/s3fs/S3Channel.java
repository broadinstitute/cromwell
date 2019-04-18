package org.lerch.s3fs;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public interface S3Channel {
    /**
     * Create the full file hierarchy for an S3 file within a temporary directory,
     * e.g. "/tmp/workflowA/workflowB/taskA/file.txt"
     * @param path S3 path to create a file for
     * @return Path on the local filesystem to the temporary file that was created
     */
    default Path createTempFile(S3Path path) throws IOException {
        return Files.createTempFile(path.getFileName().toString(), "");
    }
}
