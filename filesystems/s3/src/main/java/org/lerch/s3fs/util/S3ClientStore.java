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

    public static S3ClientStore getInstance() { return instance; }

    public S3Client getClientForBucketName( String bucketName ) {
        logger.debug("obtaining client for bucket '{}'", bucketName);
        if (bucketName == null || bucketName.trim().equals("")) {
            return DEFAULT_CLIENT;
        }

        return bucketToClientMap.computeIfAbsent(bucketName, this::generateClient);
    }

    /**
     * Generate a client for the named bucket using a default client to determine the location of the named client
     * @param bucketName the named of the bucket to make the client for
     * @return an S3 client appropriate for the region of the named bucket
     */
    protected S3Client generateClient(String bucketName){
        return this.generateClient(bucketName, DEFAULT_CLIENT);
    }

    /**
     * Generate a client for the named bucket using a default client to determine the location of the named client
     * @param bucketName the named of the bucket to make the client for
     * @param locationClient the client used to determine the location of the named bucket, recommend using DEFAULT_CLIENT
     * @return an S3 client appropriate for the region of the named bucket
     */
    protected S3Client generateClient (String bucketName, S3Client locationClient) {
        logger.debug("generating client for bucket: '{}'", bucketName);
        S3Client bucketSpecificClient;
        try {
            logger.debug("determining bucket location with getBucketLocation");
            String bucketLocation = locationClient.getBucketLocation(builder -> builder.bucket(bucketName)).locationConstraintAsString();

            bucketSpecificClient = this.clientForRegion(bucketLocation);

        } catch (S3Exception e) {
            if(e.statusCode() == 403) {
                logger.info("Cannot determine location of '{}' bucket directly. Attempting to obtain bucket location with headBucket operation", bucketName);
                try {
                    final HeadBucketResponse headBucketResponse = locationClient.headBucket(builder -> builder.bucket(bucketName));
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
            logger.warn("Unable to determine the region of bucket: '{}'. Generating a client for the current region.", bucketName);
            bucketSpecificClient = S3Client.create();
        }

        return bucketSpecificClient;
    }

    private S3Client clientForRegion(String regionString){
        // It may be useful to further cache clients for regions although at some point clients for buckets may need to be
        // specialized beyond just region end points.
        Region region = regionString.equals("") ? Region.US_EAST_1 : Region.of(regionString);
        logger.debug("bucket region is: '{}'", region.id());
        return S3Client.builder().region(region).build();
    }

}
