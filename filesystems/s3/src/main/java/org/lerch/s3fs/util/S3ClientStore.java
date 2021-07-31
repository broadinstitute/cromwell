package org.lerch.s3fs.util;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import software.amazon.awssdk.regions.Region;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.model.HeadBucketResponse;
import software.amazon.awssdk.services.s3.model.S3Exception;

import java.net.URI;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * A Singleton cache of clients for buckets configured for the region of those buckets
 */
public class S3ClientStore {

    private static final S3ClientStore instance = new S3ClientStore();

    public static final S3Client DEFAULT_CLIENT = S3Client.builder().endpointOverride(URI.create("https://s3.us-east-1.amazonaws.com")).region(Region.US_EAST_1).build();

    private final Map<String, S3Client> bucketToClientMap = Collections.synchronizedMap(new HashMap<>());

    Logger logger = LoggerFactory.getLogger("S3ClientStore");


    private S3ClientStore(){}

    public static S3ClientStore getInstance(){
        return instance;
    }

    public S3Client getClientForBucketName( String bucketName ) {
        logger.info("obtaining client for bucket '{}'", bucketName);
        if (bucketName == null || bucketName.trim().equals("")) {
            return DEFAULT_CLIENT;
        }

        return bucketToClientMap.computeIfAbsent(bucketName, this::generateClient);
    }

    private S3Client generateClient (String name) {
        logger.info("generating client for bucket: '{}'", name);
        S3Client bucketSpecificClient;
        try {
            logger.info("determining bucket location with getBucketLocation");
            String bucketLocation = DEFAULT_CLIENT.getBucketLocation(builder -> builder.bucket(name)).locationConstraintAsString();

            bucketSpecificClient = this.clientForRegion(bucketLocation);

        } catch (S3Exception e) {
            if(e.statusCode() == 403) {
                logger.info("Cannot determine bucket location directly. Attempting to obtain bucket location with headBucket operation");
                try {
                    final HeadBucketResponse headBucketResponse = DEFAULT_CLIENT.headBucket(builder -> builder.bucket(name));
                    bucketSpecificClient = this.clientForRegion(headBucketResponse.sdkHttpResponse().firstMatchingHeader("x-amz-bucket-region").orElseThrow());
                } catch (S3Exception e2) {
                    if (e2.statusCode() == 301) {
                        bucketSpecificClient = this.clientForRegion(e2.awsErrorDetails().sdkHttpResponse().firstMatchingHeader("x-amz-bucket-region").orElseThrow());
                    } else {
                        throw e2;
                    }
                }
            } else {
                throw e;
            }
        }

        if (bucketSpecificClient == null) {
            logger.warn("Unable to determine the region of bucket: '{}'", name);
            logger.warn("Generating a client for the current region");
            bucketSpecificClient = S3Client.create();
        }

        return bucketSpecificClient;
    }

    private S3Client clientForRegion(String regionString){
        // It may be useful to further cache clients for regions although at some point clients for buckets may need to be
        // specialized beyond just region end points.
        Region region = regionString.equals("") ? Region.US_EAST_1 : Region.of(regionString);
        logger.info("bucket region is: '{}'", region.id());
        return S3Client.builder().region(region).build();
    }
}
