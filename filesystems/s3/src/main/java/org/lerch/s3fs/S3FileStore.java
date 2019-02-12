package org.lerch.s3fs;

import java.io.IOException;
import java.nio.file.FileStore;
import java.nio.file.attribute.FileAttributeView;
import java.nio.file.attribute.FileStoreAttributeView;
import java.util.Date;

import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.model.Bucket;
import software.amazon.awssdk.services.s3.model.GetBucketAclRequest;
import software.amazon.awssdk.services.s3.model.HeadBucketRequest;
import software.amazon.awssdk.services.s3.model.ListBucketsRequest;
import software.amazon.awssdk.services.s3.model.NoSuchBucketException;
import software.amazon.awssdk.services.s3.model.Owner;
import com.google.common.collect.ImmutableList;

public class S3FileStore extends FileStore implements Comparable<S3FileStore> {

    private S3FileSystem fileSystem;
    private String name;

    public S3FileStore(S3FileSystem s3FileSystem, String name) {
        this.fileSystem = s3FileSystem;
        this.name = name;
    }

    @Override
    public String name() {
        return name;
    }

    @Override
    public String type() {
        return "S3Bucket";
    }

    @Override
    public boolean isReadOnly() {
        return false;
    }

    @Override
    public long getTotalSpace() throws IOException {
        return Long.MAX_VALUE;
    }

    @Override
    public long getUsableSpace() throws IOException {
        return Long.MAX_VALUE;
    }

    @Override
    public long getUnallocatedSpace() throws IOException {
        return Long.MAX_VALUE;
    }

    @Override
    public boolean supportsFileAttributeView(Class<? extends FileAttributeView> type) {
        return false;
    }

    @Override
    public boolean supportsFileAttributeView(String attributeViewName) {
        return false;
    }

    @SuppressWarnings("unchecked")
    @Override
    public <V extends FileStoreAttributeView> V getFileStoreAttributeView(Class<V> type) {
        if (type != S3FileStoreAttributeView.class)
            throw new IllegalArgumentException("FileStoreAttributeView of type '" + type.getName() + "' is not supported.");
        // TODO: bucket can come back as null here...
        Bucket buck = getBucket(name);
        Owner owner = getClient().getBucketAcl(GetBucketAclRequest.builder().bucket(name).build()).owner();
        return (V) new S3FileStoreAttributeView(Date.from(buck.creationDate()), buck.name(), owner.id(), owner.displayName());
    }

    @Override
    public Object getAttribute(String attribute) throws IOException {
        return getFileStoreAttributeView(S3FileStoreAttributeView.class).getAttribute(attribute);
    }

    public S3FileSystem getFileSystem() {
        return fileSystem;
    }

    public Bucket getBucket() {
        return getBucket(name);
    }

    private Bucket getBucket(String bucketName) {
        for (Bucket buck : getClient().listBuckets().buckets())
            if (buck.name().equals(bucketName))
                return buck;
        return null;
    }

    private boolean hasBucket(String bucketName) {
        // Originally getBucket was being used to determine presence of a bucket
        //
        // This is incorrect for two reasons:
        // 1. The list bucket operation provides buckets for which you are the owner
        //    It would not, therefore, allow you to work with buckets for which you
        //    have access but are not the owner.
        // 2. The way this information is being used later is to determine the
        //    bucket owner, which by definition, is now "you".
        // https://docs.aws.amazon.com/AmazonS3/latest/API/RESTServiceGET.html
        //
        // However, note that the revised code below now has a different permissions
        // model as HeadBucket is now required
        boolean bucket = false;
        try {
          getClient().headBucket(HeadBucketRequest.builder().bucket(bucketName).build());
          bucket = true;
        }catch(NoSuchBucketException nsbe) {}
        return bucket;
    }

    public S3Path getRootDirectory() {
        return new S3Path(fileSystem, "/" + this.name());
    }

    private S3Client getClient() {
        return fileSystem.getClient();
    }

    public Owner getOwner() {
        if (hasBucket(name))
            return getClient().getBucketAcl(GetBucketAclRequest.builder().bucket(name).build()).owner();
        // SDK v1 getS3AccountOwner uses the list buckets call, then extracts
        // the owner field (see: https://github.com/aws/aws-sdk-java/blob/4734de6fb0f80fe5768a6587aad3b9d0eaec388f/aws-java-sdk-s3/src/main/java/com/amazonaws/services/s3/model/transform/Unmarshallers.java#L48
        // and https://github.com/aws/aws-sdk-java/blob/2d15a603a96f98076f5458db49d659f296eab313/aws-java-sdk-s3/src/main/java/com/amazonaws/services/s3/AmazonS3Client.java#L926
        //
        // SDK v2 does not have that, as the SDK is mostly auto-generated based on the model files from the service, so much less custom code and helpers
        // More transparency, but we have to unwind this manually. So, here we go...
        return getClient().listBuckets(ListBucketsRequest.builder().build()).owner();
    }

    @Override
    public int compareTo(S3FileStore o) {
        if (this == o)
            return 0;
        return o.name().compareTo(name);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((fileSystem == null) ? 0 : fileSystem.hashCode());
        result = prime * result + ((name == null) ? 0 : name.hashCode());
        return result;
    }


    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (!(obj instanceof S3FileStore))
            return false;
        S3FileStore other = (S3FileStore) obj;

        if (fileSystem == null) {
            if (other.fileSystem != null)
                return false;
        } else if (!fileSystem.equals(other.fileSystem))
            return false;
        if (name == null) {
            if (other.name != null)
                return false;
        } else if (!name.equals(other.name))
            return false;
        return true;
    }
}
