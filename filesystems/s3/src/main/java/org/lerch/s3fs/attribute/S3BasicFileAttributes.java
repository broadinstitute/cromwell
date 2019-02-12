package org.lerch.s3fs.attribute;

import static java.lang.String.format;

import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileTime;

public class S3BasicFileAttributes implements BasicFileAttributes {
    private final FileTime lastModifiedTime;
    private final long size;
    private final boolean directory;
    private final boolean regularFile;
    private final String key;
    private long cacheCreated;

    public S3BasicFileAttributes(String key, FileTime lastModifiedTime, long size, boolean isDirectory, boolean isRegularFile) {
        this.key = key;
        this.lastModifiedTime = lastModifiedTime;
        this.size = size;
        this.directory = isDirectory;
        this.regularFile = isRegularFile;

        this.cacheCreated = System.currentTimeMillis();
    }

    @Override
    public FileTime lastModifiedTime() {
        return lastModifiedTime;
    }

    @Override
    public FileTime lastAccessTime() {
        return lastModifiedTime;
    }

    @Override
    public FileTime creationTime() {
        return lastModifiedTime;
    }

    @Override
    public boolean isRegularFile() {
        return regularFile;
    }

    @Override
    public boolean isDirectory() {
        return directory;
    }

    @Override
    public boolean isSymbolicLink() {
        return false;
    }

    @Override
    public boolean isOther() {
        return false;
    }

    @Override
    public long size() {
        return size;
    }

    @Override
    public Object fileKey() {
        return key;
    }

    @Override
    public String toString() {
        return format("[%s: lastModified=%s, size=%s, isDirectory=%s, isRegularFile=%s]", key, lastModifiedTime, size, directory, regularFile);
    }

    public long getCacheCreated() {
        return cacheCreated;
    }

    // for testing

    public void setCacheCreated(long time) {
        this.cacheCreated = time;
    }
}
