package org.lerch.s3fs;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

import software.amazon.awssdk.services.s3.model.CommonPrefix;
import software.amazon.awssdk.services.s3.model.ListObjectsRequest;
import software.amazon.awssdk.services.s3.model.ListObjectsResponse;
import software.amazon.awssdk.services.s3.model.S3Object;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.lerch.s3fs.util.S3Utils;

/**
 * S3 iterator over folders at first level.
 * Future versions of this class should be return the elements
 * in a incremental way when the #next() method is called.
 */
public class S3Iterator implements Iterator<Path> {
    private S3FileSystem fileSystem;
    private S3FileStore fileStore;
    private String key;
    private List<S3Path> items = Lists.newArrayList();
    private Set<S3Path> addedVirtualDirectories = Sets.newHashSet();
    private ListObjectsResponse current;
    private int cursor; // index of next element to return
    private int size;
    private boolean incremental;

    private S3Utils s3Utils = new S3Utils();

    public S3Iterator(S3Path path) {
        this(path, false);
    }

    public S3Iterator(S3Path path, boolean incremental) {
        this(path.getFileStore(), path.getKey() + (!incremental && !path.getKey().isEmpty() && !path.getKey().endsWith("/") ? "/" : ""), incremental);
    }

    public S3Iterator(S3FileStore fileStore, String key, boolean incremental) {
        // TODO: Convert to ListObjectsV2Request (marker doesn't exist in that api)
        ListObjectsRequest listObjectsRequest = buildRequest(fileStore.name(), key, incremental);

        this.fileStore = fileStore;
        this.fileSystem = fileStore.getFileSystem();
        this.key = key;
        this.current = fileSystem.getClient().listObjects(listObjectsRequest);
        this.incremental = incremental;
        loadObjects();
    }

    @Override
    public boolean hasNext() {
        return cursor != size || current.isTruncated();
    }

    @Override
    public S3Path next() {
        if (cursor == size && current.isTruncated()) {
            ListObjectsRequest request = ListObjectsRequest.builder()
                                                   .bucket(fileStore.name())
                                                   .prefix(key)
                                                   .marker(current.nextMarker())
                                                   .build();

            this.current = fileSystem.getClient().listObjects(request);
            loadObjects();
        }
        if (cursor == size)
            throw new NoSuchElementException();
        return items.get(cursor++);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    private void loadObjects() {
        this.items.clear();
        if (incremental)
            parseObjects();
        else
            parseObjectListing(key, items, current);
        this.size = items.size();
        this.cursor = 0;
    }

    private void parseObjects() {
        for (final S3Object objectSummary : current.contents()) {
            final String objectSummaryKey = objectSummary.key();
            String[] keyParts = fileSystem.key2Parts(objectSummaryKey);
            addParentPaths(keyParts);
            S3Path path = new S3Path(fileSystem, "/" + fileStore.name(), keyParts);
            if (!items.contains(path)) {
                items.add(path);
            }
        }
    }

    private void addParentPaths(String[] keyParts) {
        if (keyParts.length <= 1)
            return;
        String[] subParts = Arrays.copyOf(keyParts, keyParts.length - 1);
        List<S3Path> parentPaths = new ArrayList<>();
        while (subParts.length > 0) {
            S3Path path = new S3Path(fileSystem,  "/" + fileStore.name(), subParts);
            String prefix = current.prefix();

            String parentKey = path.getKey();
            if (prefix.length() > parentKey.length() && prefix.contains(parentKey))
                break;
            if (items.contains(path) || addedVirtualDirectories.contains(path)) {
                subParts = Arrays.copyOf(subParts, subParts.length - 1);
                continue;
            }
            parentPaths.add(path);
            addedVirtualDirectories.add(path);
            subParts = Arrays.copyOf(subParts, subParts.length - 1);
        }
        Collections.reverse(parentPaths);
        items.addAll(parentPaths);
    }


    /**
     * add to the listPath the elements at the same level that s3Path
     *
     * @param key      the uri to parse
     * @param listPath List not null list to add
     * @param current  List<S3Object> to walk
     */
    private void parseObjectListing(String key, List<S3Path> listPath, ListObjectsResponse current) {
        for (CommonPrefix commonPrefix : current.commonPrefixes()) {
            if (!commonPrefix.prefix().equals("/")) {
                listPath.add(new S3Path(fileSystem,  "/" + fileStore.name(), fileSystem.key2Parts(commonPrefix.prefix())));
            }
        }
        // TODO: figure our a way to efficiently preprocess commonPrefix basicFileAttributes
        for (final S3Object objectSummary : current.contents()) {
            final String objectSummaryKey = objectSummary.key();
            // we only want the first level
            String immediateDescendantKey = getImmediateDescendant(key, objectSummaryKey);
            if (immediateDescendantKey != null) {
                S3Path descendentPart = new S3Path(fileSystem,  "/" + fileStore.name(), fileSystem.key2Parts(immediateDescendantKey));
                descendentPart.setFileAttributes(s3Utils.toS3FileAttributes(objectSummary, descendentPart.getKey()));
                if (!listPath.contains(descendentPart)) {
                    listPath.add(descendentPart);
                }
            }
        }
    }

    /**
     * The current #buildRequest() get all subdirectories and her content.
     * This method filter the keyChild and check if is a inmediate
     * descendant of the keyParent parameter
     *
     * @param keyParent String
     * @param keyChild  String
     * @return String parsed
     * or null when the keyChild and keyParent are the same and not have to be returned
     */
    private String getImmediateDescendant(String keyParent, String keyChild) {
        keyParent = deleteExtraPath(keyParent);
        keyChild = deleteExtraPath(keyChild);
        final int parentLen = keyParent.length();
        final String childWithoutParent = deleteExtraPath(keyChild.substring(parentLen));
        String[] parts = childWithoutParent.split("/");
        if (parts.length > 0 && !parts[0].isEmpty())
            return keyParent + "/" + parts[0];
        return null;

    }

    private String deleteExtraPath(String keyChild) {
        if (keyChild.startsWith("/"))
            keyChild = keyChild.substring(1);
        if (keyChild.endsWith("/"))
            keyChild = keyChild.substring(0, keyChild.length() - 1);
        return keyChild;
    }


    ListObjectsRequest buildRequest(String bucketName, String key, boolean incremental) {
        return buildRequest(bucketName, key, incremental, null);
    }

    ListObjectsRequest buildRequest(String bucketName, String key, boolean incremental, Integer maxKeys) {
        ListObjectsRequest.Builder builder = ListObjectsRequest.builder();
        builder.bucket(bucketName)
               .prefix(key)
               .maxKeys(maxKeys);

        if (!incremental)
            builder.marker(key)
                   .delimiter("/");
        return builder.build();
    }
}
