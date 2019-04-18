package org.lerch.s3fs;

import static java.lang.String.format;

import java.io.ByteArrayInputStream;
import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.*;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import org.apache.tika.Tika;

import software.amazon.awssdk.core.sync.RequestBody;
import software.amazon.awssdk.core.ResponseBytes;
import software.amazon.awssdk.services.s3.S3Client;
import software.amazon.awssdk.services.s3.model.GetObjectResponse;
import software.amazon.awssdk.services.s3.model.GetObjectRequest;
import software.amazon.awssdk.services.s3.model.PutObjectRequest;
import software.amazon.awssdk.services.s3.model.S3Object;

public class S3SeekableByteChannel implements SeekableByteChannel, S3Channel {

    private S3Path path;
    private Set<? extends OpenOption> options;
    private SeekableByteChannel seekable;
    private Path tempFile;

    /**
     * Open or creates a file, returning a seekable byte channel
     *
     * @param path    the path open or create
     * @param options options specifying how the file is opened
     * @throws IOException if an I/O error occurs
     */
    public S3SeekableByteChannel(S3Path path, Set<? extends OpenOption> options) throws IOException {
        this.path = path;
        this.options = Collections.unmodifiableSet(new HashSet<>(options));
        String key = path.getKey();
        boolean exists = path.getFileSystem().provider().exists(path);

        if (exists && this.options.contains(StandardOpenOption.CREATE_NEW))
            throw new FileAlreadyExistsException(format("target already exists: %s", path));
        else if (!exists && !this.options.contains(StandardOpenOption.CREATE_NEW) &&
                !this.options.contains(StandardOpenOption.CREATE))
            throw new NoSuchFileException(format("target not exists: %s", path));

        tempFile = createTempFile(path);
        boolean removeTempFile = true;
        try {
            if (exists) {
                try (InputStream byteStream = path.getFileSystem().getClient()
                      .getObject(GetObjectRequest.builder().bucket(path.getFileStore().getBucket().name()).key(key).build())) {
                   Files.copy(byteStream, tempFile, StandardCopyOption.REPLACE_EXISTING);
                }
            }

            Set<? extends OpenOption> seekOptions = new HashSet<>(this.options);
            seekOptions.remove(StandardOpenOption.CREATE_NEW);
            seekable = Files.newByteChannel(tempFile, seekOptions);
            removeTempFile = false;
        } finally {
            if (removeTempFile) {
                Files.deleteIfExists(tempFile);
            }
        }
    }

    @Override
    public boolean isOpen() {
        return seekable.isOpen();
    }

    @Override
    public void close() throws IOException {
        try {
            if (!seekable.isOpen())
                return;

            seekable.close();

            if (options.contains(StandardOpenOption.DELETE_ON_CLOSE)) {
                path.getFileSystem().provider().delete(path);
                return;
            }

            if (options.contains(StandardOpenOption.READ) && options.size() == 1) {
                return;
            }

            sync();

        } finally {
            Files.deleteIfExists(tempFile);
        }
    }

    /**
     * try to sync the temp file with the remote s3 path.
     *
     * @throws IOException if the tempFile fails to open a newInputStream
     */
    protected void sync() throws IOException {
        try (InputStream stream = new BufferedInputStream(Files.newInputStream(tempFile))) {
            PutObjectRequest.Builder builder = PutObjectRequest.builder();
            long length = Files.size(tempFile);
            builder.contentLength(length);
            if (path.getFileName() != null) {
                builder.contentType(new Tika().detect(stream, path.getFileName().toString()));
            }

            builder.bucket(path.getFileStore().name());
            builder.key(path.getKey());

            S3Client client = path.getFileSystem().getClient();

            client.putObject(builder.build(), RequestBody.fromInputStream(stream, length));
        }
    }

    @Override
    public int write(ByteBuffer src) throws IOException {
        return seekable.write(src);
    }

    @Override
    public SeekableByteChannel truncate(long size) throws IOException {
        return seekable.truncate(size);
    }

    @Override
    public long size() throws IOException {
        return seekable.size();
    }

    @Override
    public int read(ByteBuffer dst) throws IOException {
        return seekable.read(dst);
    }

    @Override
    public SeekableByteChannel position(long newPosition) throws IOException {
        return seekable.position(newPosition);
    }

    @Override
    public long position() throws IOException {
        return seekable.position();
    }
}
