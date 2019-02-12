package org.lerch.s3fs.attribute;

import org.lerch.s3fs.S3Path;

import java.io.IOException;
import java.nio.file.attribute.*;
import java.util.Set;


public class S3PosixFileAttributeView implements PosixFileAttributeView {

    private S3Path s3Path;
    private PosixFileAttributes posixFileAttributes;

    public S3PosixFileAttributeView(S3Path s3Path) {
        this.s3Path = s3Path;
    }

    @Override
    public String name() {
        return "posix";
    }

    @Override
    public PosixFileAttributes readAttributes() throws IOException {
        return read();
    }

    @Override
    public UserPrincipal getOwner() throws IOException {
        return read().owner();
    }

    @Override
    public void setOwner(UserPrincipal owner) throws IOException {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setPermissions(Set<PosixFilePermission> perms) throws IOException {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setGroup(GroupPrincipal group) throws IOException {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setTimes(FileTime lastModifiedTime, FileTime lastAccessTime, FileTime createTime) throws IOException {
        // not implemented
    }

    public PosixFileAttributes read() throws IOException {
        if (posixFileAttributes == null) {
            posixFileAttributes = s3Path.getFileSystem()
                    .provider().readAttributes(s3Path, PosixFileAttributes.class);
        }
        return posixFileAttributes;
    }
}
