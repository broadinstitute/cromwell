package org.lerch.s3fs.util;

import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.model.HeadObjectRequest;
import software.amazon.awssdk.services.s3.model.HeadObjectResponse;
import software.amazon.awssdk.services.s3.model.GetObjectAclRequest;
import software.amazon.awssdk.services.s3.model.GetObjectAclResponse;
import software.amazon.awssdk.services.s3.model.Grant;
import software.amazon.awssdk.services.s3.model.ListObjectsV2Request;
import software.amazon.awssdk.services.s3.model.ListObjectsV2Response;
import software.amazon.awssdk.services.s3.model.Owner;
import software.amazon.awssdk.services.s3.model.Permission;
import software.amazon.awssdk.services.s3.model.S3Object;
import software.amazon.awssdk.services.s3.model.S3Exception;
import com.google.common.collect.Sets;
import org.lerch.s3fs.attribute.S3BasicFileAttributes;
import org.lerch.s3fs.S3Path;
import org.lerch.s3fs.attribute.S3PosixFileAttributes;
import org.lerch.s3fs.attribute.S3UserPrincipal;

import java.nio.file.NoSuchFileException;
import java.nio.file.attribute.FileTime;
import java.nio.file.attribute.PosixFilePermission;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.TimeUnit;

/**
 * Utilities to work with Amazon S3 Objects.
 */
public class S3Utils {

    /**
     * Get the {@link S3Object} that represent this Path or her first child if this path not exists
     *
     * @param s3Path {@link S3Path}
     * @return {@link S3Object}
     * @throws NoSuchFileException if not found the path and any child
     */
    public S3Object getS3ObjectSummary(S3Path s3Path) throws NoSuchFileException {
        String key = s3Path.getKey();
        String bucketName = s3Path.getFileStore().name();
        S3Client client = s3Path.getFileSystem().getClient();
        // try to find the element with the current key (maybe with end slash or maybe not.)
        try {
            HeadObjectResponse metadata = client.headObject(HeadObjectRequest.builder().bucket(bucketName).key(key).build());
            GetObjectAclResponse acl = client.getObjectAcl(GetObjectAclRequest.builder().bucket(bucketName).key(key).build());
            S3Object.Builder builder = S3Object.builder();

            builder
                .key(key)
                .lastModified(metadata.lastModified())
                .eTag(metadata.eTag())
                .owner(acl.owner())
                .size(metadata.contentLength())
                .storageClass(metadata.storageClassAsString());

            return builder.build();
        } catch (S3Exception e) {
            if (e.statusCode() != 404)
                throw e;
        }

        // if not found (404 err) with the original key.
        // try to find the elment as a directory.
        try {
            // is a virtual directory
            ListObjectsV2Request.Builder request = ListObjectsV2Request.builder();
            request.bucket(bucketName);
            String keyFolder = key;
            if (!keyFolder.endsWith("/")) {
                keyFolder += "/";
            }
            request.prefix(keyFolder);
            request.maxKeys(1);
            ListObjectsV2Response current = client.listObjectsV2(request.build());
            if (!current.contents().isEmpty())
                return current.contents().get(0);
        } catch (Exception e) {
            //
        }
        throw new NoSuchFileException(bucketName + S3Path.PATH_SEPARATOR + key);
    }

    /**
     * getS3FileAttributes for the s3Path
     *
     * @param s3Path S3Path mandatory not null
     * @return S3FileAttributes never null
     */
    public S3BasicFileAttributes getS3FileAttributes(S3Path s3Path) throws NoSuchFileException {
        S3Object objectSummary = getS3ObjectSummary(s3Path);
        return toS3FileAttributes(objectSummary, s3Path.getKey());
    }

    /**
     * get the S3PosixFileAttributes for a S3Path
     * @param s3Path Path mandatory not null
     * @return S3PosixFileAttributes never null
     * @throws NoSuchFileException if the Path doesnt exists
     */
    public S3PosixFileAttributes getS3PosixFileAttributes(S3Path s3Path) throws NoSuchFileException {
        S3Object objectSummary = getS3ObjectSummary(s3Path);

        String key = s3Path.getKey();
        String bucketName = s3Path.getFileStore().name();

        S3BasicFileAttributes attrs = toS3FileAttributes(objectSummary, key);
        S3UserPrincipal userPrincipal = null;
        Set<PosixFilePermission> permissions = null;

        if (!attrs.isDirectory()) {
            S3Client client = s3Path.getFileSystem().getClient();
            GetObjectAclResponse acl = client.getObjectAcl(GetObjectAclRequest.builder().bucket(bucketName).key(key).build());
            Owner owner = acl.owner();

            userPrincipal = new S3UserPrincipal(owner.id() + ":" + owner.displayName());
            permissions = toPosixFilePermissions(acl.grants());
        }

        return new S3PosixFileAttributes((String)attrs.fileKey(), attrs.lastModifiedTime(),
                attrs.size(), attrs.isDirectory(), attrs.isRegularFile(), userPrincipal, null, permissions);
    }


    /**
     * transform software.amazon.awssdk.services.s3.model.Grant to java.nio.file.attribute.PosixFilePermission
     * @see #toPosixFilePermission(Permission)
     * @param grants Set grants mandatory, must be not null
     * @return Set PosixFilePermission never null
     */
    public Set<PosixFilePermission> toPosixFilePermissions(List<Grant> grants) {
        Set<PosixFilePermission> filePermissions = new HashSet<>();
        for (Grant grant : grants) {
            filePermissions.addAll(toPosixFilePermission(grant.permission()));
        }

        return filePermissions;
    }

    /**
     * transform a software.amazon.awssdk.services.s3.model.Permission to a java.nio.file.attribute.PosixFilePermission
     * We use the follow rules:
     * - transform only to the Owner permission, S3 doesnt have concepts like owner, group or other so we map only to owner.
     * - ACP is a special permission: WriteAcp are mapped to Owner execute permission and ReadAcp are mapped to owner read
     * @param permission Permission to map, mandatory must be not null
     * @return Set PosixFilePermission never null
     */
    public Set<PosixFilePermission> toPosixFilePermission(Permission permission){
        switch (permission) {
            case FULL_CONTROL:
                return Sets.newHashSet(PosixFilePermission.OWNER_EXECUTE,
                        PosixFilePermission.OWNER_READ,
                        PosixFilePermission.OWNER_WRITE);
            case WRITE:
                return Sets.newHashSet(PosixFilePermission.OWNER_WRITE);
            case READ:
                return Sets.newHashSet(PosixFilePermission.OWNER_READ);
            case READ_ACP:
                return Sets.newHashSet(PosixFilePermission.OWNER_READ);
            case WRITE_ACP:
                return Sets.newHashSet(PosixFilePermission.OWNER_EXECUTE);
        }
        throw new IllegalStateException("Unknown Permission: " + permission);
    }

    /**
     * transform S3ObjectSummary to S3FileAttributes
     *
     * @param objectSummary S3ObjectSummary mandatory not null, the real objectSummary with
     *                      exactly the same key than the key param or the immediate descendant
     *                      if it is a virtual directory
     * @param key           String the real key that can be exactly equal than the objectSummary or
     * @return S3FileAttributes
     */
    public S3BasicFileAttributes toS3FileAttributes(S3Object objectSummary, String key) {
        // parse the data to BasicFileAttributes.
        FileTime lastModifiedTime = null;
        if (objectSummary.lastModified() != null) {
            lastModifiedTime = FileTime.from(objectSummary.lastModified());
        }
        long size = objectSummary.size();
        boolean directory = false;
        boolean regularFile = false;
        String resolvedKey = objectSummary.key();
        // check if is a directory and exists the key of this directory at amazon s3
        if (key.endsWith("/") && resolvedKey.equals(key) ||
                resolvedKey.equals(key + "/")) {
            directory = true;
        } else if (key.isEmpty()) { // is a bucket (no key)
            directory = true;
            resolvedKey = "/";
        } else if (!resolvedKey.equals(key) && resolvedKey.startsWith(key)) { // is a directory but not exists at amazon s3
            directory = true;
            // no metadata, we fake one
            size = 0;
            // delete extra part
            resolvedKey = key + "/";
        } else {
            regularFile = true;
        }
        return new S3BasicFileAttributes(resolvedKey, lastModifiedTime, size, directory, regularFile);
    }
}
