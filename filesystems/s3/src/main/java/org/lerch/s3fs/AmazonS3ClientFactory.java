package org.lerch.s3fs;

import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.S3ClientBuilder;

public class AmazonS3ClientFactory extends AmazonS3Factory {

    @Override
    protected S3Client createS3Client(S3ClientBuilder builder) {
        return builder.build();
    }
}
