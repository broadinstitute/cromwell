package org.lerch.s3fs;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public interface S3Channel {
    default Path createTempFile(S3Path path) throws IOException {
        return Files.createTempFile(path.getFileName().toString(), "");
    }
}
