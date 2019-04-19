package org.lerch.s3fs.attribute;

import java.nio.file.attribute.*;
import java.util.Set;

import static java.lang.String.format;

public class S3PosixFileAttributes extends S3BasicFileAttributes implements PosixFileAttributes  {

    private UserPrincipal userPrincipal;
    private GroupPrincipal groupPrincipal;
    private Set<PosixFilePermission> posixFilePermissions;

    public S3PosixFileAttributes(String key, FileTime lastModifiedTime, long size, boolean isDirectory, boolean isRegularFile, UserPrincipal userPrincipal, GroupPrincipal groupPrincipal, Set<PosixFilePermission> posixFilePermissionSet) {

        super(key, lastModifiedTime, size, isDirectory, isRegularFile);

        this.userPrincipal = userPrincipal;
        this.groupPrincipal = groupPrincipal;
        this.posixFilePermissions = posixFilePermissionSet;
    }

    @Override
    public UserPrincipal owner() {
        return this.userPrincipal;
    }

    @Override
    public GroupPrincipal group() {
        return this.groupPrincipal;
    }

    @Override
    public Set<PosixFilePermission> permissions() {
        return this.posixFilePermissions;
    }
}
