package org.lerch.s3fs.util;

import com.amazonaws.auth.AWSStaticCredentialsProvider;
import com.amazonaws.auth.BasicAWSCredentials;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3ClientBuilder;
import software.amazon.awssdk.auth.credentials.AwsCredentials;
import software.amazon.awssdk.regions.Region;

import java.util.Objects;

public class AmazonS3ClientProvider {
    private final AwsCredentials credentials;
    private final Region region;
    private static AmazonS3ClientProvider instance;

    public static void init(AwsCredentials credentials, Region region) {
        if (Objects.isNull(instance))
            instance = new AmazonS3ClientProvider(credentials, region);
    }

    public static AmazonS3 buildAmazonS3Client() throws AssertionError {
        if (Objects.isNull(instance))
            throw new AssertionError("AmazonS3ClientProvider instance should no be null!");

        final AwsCredentials creds = instance.credentials;
        final Region reg = instance.region;

        final BasicAWSCredentials basicAWSCredentials = new BasicAWSCredentials(creds.accessKeyId(),
                creds.secretAccessKey());
        return AmazonS3ClientBuilder.standard()
                .withCredentials(new AWSStaticCredentialsProvider(basicAWSCredentials))
                .withRegion(reg.id())
                .build();
    }

    private AmazonS3ClientProvider(AwsCredentials credentials, Region region) {
        this.credentials = credentials;
        this.region = region;
    }
}
