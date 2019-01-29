package org.lerch.s3fs;

import static java.lang.String.format;

import java.nio.file.AccessDeniedException;
import java.nio.file.AccessMode;
import java.util.EnumSet;

import software.amazon.awssdk.services.s3.model.Grant;
import software.amazon.awssdk.services.s3.model.Owner;
import software.amazon.awssdk.services.s3.model.Permission;

public class S3AccessControlList {
    private String fileStoreName;
    private String key;
    private Iterable<Grant> grants;
    private Owner owner;

    public S3AccessControlList(String fileStoreName, String key, Iterable<Grant> grants, Owner owner) {
        this.fileStoreName = fileStoreName;
        this.grants = grants;
        this.key = key;
        this.owner = owner;
    }

    public String getKey() {
        return key;
    }

    /**
     * have almost one of the permission set in the parameter permissions
     *
     * @param permissions almost one
     * @return
     */
    private boolean hasPermission(EnumSet<Permission> permissions) {
        for (Grant grant : grants)
            if (grant.grantee().id().equals(owner.id()) && permissions.contains(grant.permission()))
                return true;
        return false;
    }

    public void checkAccess(AccessMode[] modes) throws AccessDeniedException {
        for (AccessMode accessMode : modes) {
            switch (accessMode) {
                case EXECUTE:
                    throw new AccessDeniedException(fileName(), null, "file is not executable");
                case READ:
                    if (!hasPermission(EnumSet.of(Permission.FULL_CONTROL, Permission.READ)))
                        throw new AccessDeniedException(fileName(), null, "file is not readable");
                    break;
                case WRITE:
                    if (!hasPermission(EnumSet.of(Permission.FULL_CONTROL, Permission.WRITE)))
                        throw new AccessDeniedException(fileName(), null, format("bucket '%s' is not writable", fileStoreName));
                    break;
            }
        }
    }

    private String fileName() {
        return fileStoreName + S3Path.PATH_SEPARATOR + key;
    }
}
