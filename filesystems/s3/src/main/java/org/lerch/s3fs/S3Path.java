package org.lerch.s3fs;

import com.google.common.base.*;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import org.lerch.s3fs.attribute.S3BasicFileAttributes;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URL;
import java.net.URLDecoder;
import java.nio.file.*;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import static com.google.common.collect.Iterables.*;
import static java.lang.String.format;

public class S3Path implements Path {

    public static final String PATH_SEPARATOR = "/";
    /**
     * S3FileStore which represents the Bucket this path resides in.
     */
    private final S3FileStore fileStore;

    /**
     * URI not encoded
     * Is the key for AmazonS3
     */
    private String uri;

    /**
     * actual filesystem
     */
    private S3FileSystem fileSystem;

    /**
     * S3BasicFileAttributes cache
     */
    private S3BasicFileAttributes fileAttributes;

    /**
     * Build an S3Path from path segments. '/' are stripped from each segment.
     *
     * @param fileSystem S3FileSystem
     * @param first should be start with a '/' and is the bucket name
     * @param more  directories and files
     */
    public S3Path(S3FileSystem fileSystem, String first, String... more) {

        Preconditions.checkArgument(first != null, "first path must be not null");
        Preconditions.checkArgument(!first.startsWith("//"), "first path doesnt start with '//'. Miss bucket");
        // see tests org.lerch.s3fs.Path.EndsWithTest#endsWithRelativeBlankAbsolute()
        // Preconditions.checkArgument(!first.isEmpty(), "first path must be not empty");

        boolean hasBucket = first.startsWith("/");


        List<String> pathsURI = Lists
                .newArrayList(Splitter.on(PATH_SEPARATOR)
                        .omitEmptyStrings()
                        .split(first));

        if (hasBucket) { // absolute path

            Preconditions.checkArgument(pathsURI.size() >= 1, "path must start with bucket name");
            Preconditions.checkArgument(!pathsURI.get(0).isEmpty(), "bucket name must be not empty");
            String bucket = pathsURI.get(0);
            this.fileStore = new S3FileStore(fileSystem, bucket);
            // the filestore is not part of the uri
            pathsURI.remove(0);
        }
        else {
            // relative uri
            this.fileStore = null;
        }

        StringBuilder uriBuilder = new StringBuilder();
        if (hasBucket) {
            uriBuilder.append(PATH_SEPARATOR);
        }
        for (String path : pathsURI) {
            uriBuilder.append(path + PATH_SEPARATOR);
        }
        if (more != null) {
            for (String path : more) {
                uriBuilder.append(path + PATH_SEPARATOR);
            }
        }
        this.uri = normalizeURI(uriBuilder.toString());
        // remove last PATH_SEPARATOR
        if (!first.isEmpty() &&
                // only first param and not ended with PATH_SEPARATOR
                ((!first.endsWith(PATH_SEPARATOR) && (more == null || more.length == 0))
                // we have more param and not ended with PATH_SEPARATOR
                || more != null &&  more.length > 0 && !more[more.length-1].endsWith(PATH_SEPARATOR))) {
            this.uri = this.uri.substring(0, this.uri.length() - 1);
        }

        this.fileSystem = fileSystem;
    }

    /**
     * Remove duplicated slash
     */
    private String normalizeURI(String uri) {
        return uri.replace("//", "/");
    }

    public S3FileStore getFileStore() {
        return fileStore;
    }

    /**
     * key for amazon without final slash.
     * <b>note:</b> the final slash need to be added to save a directory (Amazon s3 spec)
     *
     * @return the key for AmazonS3Client
     */
    public String getKey() {

        String key = this.uri;

        if (key.startsWith("/")) {
            key = key.substring(1, key.length());
        }

        return key;
    }

    @Override
    public S3FileSystem getFileSystem() {
        return this.fileSystem;
    }

    @Override
    public boolean isAbsolute() {
        return fileStore != null;
    }

    @Override
    public Path getRoot() {
        if (isAbsolute()) {
            return new S3Path(fileSystem, PATH_SEPARATOR + fileStore.name() + PATH_SEPARATOR);
        }

        return null;
    }

    @Override
    public Path getFileName() {
        List<String> paths = uriToList();
        if (paths.isEmpty()) {
            // get FileName of root directory is null
            return null;
        }
        String filename = paths.get(paths.size()-1);
        return new S3Path(fileSystem, filename);
    }

    @Override
    public Path getParent() {
        // bucket is not present in the parts
        if (uri.isEmpty()) {
            return null;
        }

        String newUri = this.uri;

        if (this.uri.endsWith("/")) {
            newUri = this.uri.substring(0, this.uri.length()-1);
        }
        int lastPathSeparatorPosition = newUri.lastIndexOf(PATH_SEPARATOR);

        if (lastPathSeparatorPosition == -1) {
            return null;
        }

        newUri = uri.substring(0, lastPathSeparatorPosition + 1);

        if (newUri.isEmpty())
            return null;

        String filestore = isAbsolute() ? PATH_SEPARATOR + fileStore.name() + PATH_SEPARATOR : "";

        return new S3Path(fileSystem, filestore + newUri);
    }

    @Override
    public int getNameCount() {
        return uriToList().size();
    }

    @Override
    public Path getName(int index) {

        List<String> paths = uriToList();

        if (index < 0 || index >= paths.size()) {
            throw new IllegalArgumentException("index out of range");
        }

        String path = paths.get(index);
        StringBuilder pathsBuilder = new StringBuilder();
        if (isAbsolute() && index == 0) {
            pathsBuilder.append(PATH_SEPARATOR + fileStore.name() + PATH_SEPARATOR);
        }
        pathsBuilder.append(path);

        if (index < paths.size() - 1) {
            pathsBuilder.append(PATH_SEPARATOR);
        }

        // if is the last path, check if end with path separator
        if (index == paths.size() - 1 && this.uri.endsWith(PATH_SEPARATOR)) {
            pathsBuilder.append(PATH_SEPARATOR);
        }

        return new S3Path(fileSystem, pathsBuilder.toString());
    }

    private List<String> uriToList() {
        return Splitter.on(PATH_SEPARATOR).omitEmptyStrings().splitToList(this.uri);
    }

    /**
     * The bucket name not count
     */
    @Override
    public Path subpath(int beginIndex, int endIndex) {

        List<String> paths = uriToList();

        if (beginIndex < 0 || endIndex > paths.size()) {
            throw new IllegalArgumentException("index out of range");
        }

        List<String> pathSubList = paths.subList(beginIndex, endIndex);
        StringBuilder pathsStringBuilder = new StringBuilder();

        // build path string

        if (this.isAbsolute() && beginIndex == 0) {
            pathsStringBuilder.append(PATH_SEPARATOR + fileStore.name() + PATH_SEPARATOR);
        }
        for (String path : pathSubList) {
            pathsStringBuilder.append(path).append(PATH_SEPARATOR);
        }
        String pathsResult = pathsStringBuilder.toString();
        // if the uri doesnt have last PATH_SEPARATOR we must remove it.
        if (endIndex == paths.size() && !this.uri.endsWith(PATH_SEPARATOR)) {
            pathsResult = pathsResult.substring(0, pathsResult.length() - 1);
        }

        return new S3Path(fileSystem, pathsResult);
    }

    @Override
    public boolean startsWith(Path other) {

        if (other.getNameCount() > this.getNameCount()) {
            return false;
        }

        if (!(other instanceof S3Path)) {
            return false;
        }

        if (this.isAbsolute() && !other.isAbsolute()) {
            return false;
        }

        S3Path path = (S3Path) other;

        if (this.isAbsolute() && other.isAbsolute() &&
                !this.fileStore.name().equals(path.fileStore.name())) {
            return false;
        }

        if (path.uri.isEmpty() && !this.uri.isEmpty()) {
            return false;
        }

        List<String> pathsOther = path.uriToList();
        List<String> paths = this.uriToList();


        for (int i = 0; i < pathsOther.size(); i++) {
            if (!pathsOther.get(i).equals(paths.get(i))) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean startsWith(String path) {
        S3Path other = new S3Path(this.fileSystem, path);
        return this.startsWith(other);
    }

    @Override
    public boolean endsWith(Path other) {
        if (other.getNameCount() > this.getNameCount()) {
            return false;
        }
        // empty
        if (other.getNameCount() == 0 && this.getNameCount() != 0) {
            return false;
        }

        if (!(other instanceof S3Path)) {
            return false;
        }

        S3Path path = (S3Path) other;

        if ((path.getFileStore() != null && !path.getFileStore().equals(this.getFileStore())) || (path.getFileStore() != null && this.getFileStore() == null)) {
            return false;
        }

        // check subkeys

        List<String> pathsOther = path.uriToList();
        List<String> paths = this.uriToList();


        int i = pathsOther.size() - 1;
        int j = paths.size() - 1;
        for (; i >= 0 && j >= 0; ) {

            if (!pathsOther.get(i).equals(paths.get(j))) {
                return false;
            }
            i--;
            j--;
        }
        return true;
    }

    @Override
    public boolean endsWith(String other) {
        return this.endsWith(new S3Path(this.fileSystem, other));
    }

    @Override
    public Path normalize() {
        return this;
    }

    @Override
    public Path resolve(Path other) {
        if (other.isAbsolute()) {
            Preconditions.checkArgument(other instanceof S3Path, "other must be an instance of %s", S3Path.class.getName());
            return other;
        }

        S3Path otherS3Path = (S3Path) other;
        StringBuilder pathBuilder = new StringBuilder();

        if (this.isAbsolute()) {
            pathBuilder.append(PATH_SEPARATOR + this.fileStore.name() + PATH_SEPARATOR);
        }
        pathBuilder.append(this.uri);
        if (!otherS3Path.uri.isEmpty())
            pathBuilder.append(PATH_SEPARATOR + otherS3Path.uri);

        return new S3Path(this.fileSystem, pathBuilder.toString());
    }

    @Override
    public Path resolve(String other) {
        return resolve(new S3Path(this.getFileSystem(), other));
    }

    @Override
    public Path resolveSibling(Path other) {
        Preconditions.checkArgument(other instanceof S3Path, "other must be an instance of %s", S3Path.class.getName());

        S3Path s3Path = (S3Path) other;

        Path parent = getParent();

        if (parent == null || s3Path.isAbsolute()) {
            return s3Path;
        }

        List<String> othersPaths = s3Path.uriToList();

        if (othersPaths.isEmpty()) { // other is relative and empty
            return parent;
        }

        List<String> paths = this.uriToList();

        StringBuilder pathBuilder = new StringBuilder();
        String lastPath = othersPaths.get(othersPaths.size() - 1);
        if (isAbsolute()) {
            pathBuilder.append(PATH_SEPARATOR + fileStore.name() + PATH_SEPARATOR);
        }
        for (String path : concat(paths.subList(0, paths.size() - 1), othersPaths)) {
            pathBuilder.append(path);
            if (!lastPath.equals(path) || s3Path.uri.endsWith(PATH_SEPARATOR)) {
                pathBuilder.append(PATH_SEPARATOR);
            }
        }

        return new S3Path(fileSystem, pathBuilder.toString());
    }

    @Override
    public Path resolveSibling(String other) {
        return resolveSibling(new S3Path(this.getFileSystem(), other));
    }

    @Override
    public Path relativize(Path other) {
        Preconditions.checkArgument(other instanceof S3Path, "other must be an instance of %s", S3Path.class.getName());
        S3Path s3Path = (S3Path) other;

        if (this.equals(other)) {
            return new S3Path(this.getFileSystem(), "");
        }

        Preconditions.checkArgument(isAbsolute(), "Path is already relative: %s", this);
        Preconditions.checkArgument(s3Path.isAbsolute(), "Cannot relativize against a relative path: %s", s3Path);
        Preconditions.checkArgument(fileStore.equals(s3Path.getFileStore()), "Cannot relativize paths with different buckets: '%s', '%s'", this, other);
        // Preconditions.checkArgument(parts.size() <= s3Path.parts.size(), "Cannot relativize against a parent path: '%s', '%s'", this, other);

        String uriPath = decode(URI.create(encode(this.uri)).relativize(URI.create(encode(s3Path.uri))));
        return new S3Path(fileSystem, uriPath);
    }

    /**
     * Examples:
     *
     * Relative:
     * --------
     * NO use fileSystem and not used fileStore.
     * - path/file
     *
     * Absolute:
     * --------
     * Use the fileSystem to get the host and the filestore to get the first path (in the future the filestore can be attached to the host)
     * http://docs.aws.amazon.com/AmazonS3/latest/dev/BucketRestrictions.html
     * - s3://AMAZONACCESSKEY@s3.amazonaws.com/bucket/path/file
     * - s3://AMAZONACCESSKEY@bucket.s3.amazonaws.com/path/file
     * - s3://s3-aws-region.amazonaws.com/bucket/path/file
     *
     * @return URI never null
     */
    @Override
    public URI toUri() {

        String uri = encode(this.uri);
        // absolute
        if (this.isAbsolute()) {
            StringBuilder builder = new StringBuilder();
            builder.append(fileSystem.getKey());
            builder.append(PATH_SEPARATOR + fileStore.name() + PATH_SEPARATOR);
            builder.append(uri);
            return URI.create("s3://" + normalizeURI(builder.toString()));
        }
        else {
            return URI.create(this.uri);
        }
    }

    @Override
    public Path toAbsolutePath() {
        if (isAbsolute()) {
            return this;
        }

        throw new IllegalStateException(format("Relative path cannot be made absolute: %s", this));
    }

    @Override
    public Path toRealPath(LinkOption... options) throws IOException {
        return toAbsolutePath();
    }

    @Override
    public File toFile() {
        throw new UnsupportedOperationException();
    }

    @Override
    public WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) throws IOException {
        throw new UnsupportedOperationException();
    }

    @Override
    public WatchKey register(WatchService watcher, WatchEvent.Kind<?>... events) throws IOException {
        throw new UnsupportedOperationException();
    }

    @Override
    public Iterator<Path> iterator() {

        ImmutableList.Builder<Path> builder = ImmutableList.builder();

        if (isAbsolute()) {
            builder.add(new S3Path(fileSystem, PATH_SEPARATOR + fileStore.name() + PATH_SEPARATOR));
        }

        List<String> paths = uriToList();

        if (uriToList().isEmpty())
            return builder.build().iterator();

        String lastPath = paths.get(paths.size() - 1);

        for (String path : uriToList()) {
            String pathFinal = path + PATH_SEPARATOR;
            if (path.equals(lastPath) && !lastPath.endsWith(PATH_SEPARATOR)) {
                pathFinal = pathFinal.substring(0, pathFinal.length() - 1);
            }
            builder.add(new S3Path(fileSystem, pathFinal));
        }

        return builder.build().iterator();
    }

    @Override
    public int compareTo(Path other) {
        return toString().compareTo(other.toString());
    }

    @Override
    public String toString() {
        return toUri().toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o)
            return true;
        if (o == null || getClass() != o.getClass())
            return false;

        S3Path path = (S3Path) o;
        if (fileStore != null ? !fileStore.equals(path.fileStore) : path.fileStore != null)
            return false;
        if (!uri.equals(path.uri))
            return false;
        return true;
    }

    @Override
    public int hashCode() {
        int result = fileStore != null ? fileStore.name().hashCode() : 0;
        result = 31 * result + uri.hashCode();
        return result;
    }

    /**
     * Encode special URI characters for path.
     * @param uri String the uri path
     * @return String
     */
    private String encode(String uri) {
        // remove special case URI starting with //
        uri = uri.replace("//", "/");
        uri = uri.replaceAll(" ", "%20");
        return uri;
    }

    /**
     * Decode uri special characters
     *
     * @param uri URI mandatory
     * @return String decoded
     */
    private String decode(URI uri) {
        try {
            return URLDecoder.decode(uri.toString(), "UTF-8");
        }
        catch (UnsupportedEncodingException e) {
            throw new IllegalStateException("Error decoding key: " + this.uri, e);
        }
    }

    public S3BasicFileAttributes getFileAttributes() {
        return fileAttributes;
    }

    public void setFileAttributes(S3BasicFileAttributes fileAttributes) {
        this.fileAttributes = fileAttributes;
    }
}
