package org.lerch.s3fs;

public class S3FileSystemConfigurationException extends RuntimeException {
    private static final long serialVersionUID = 1L;

    public S3FileSystemConfigurationException(String message, Throwable cause) {
        super(message, cause);
    }
}