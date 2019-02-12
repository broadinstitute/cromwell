package org.lerch.s3fs.util;


import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.PosixFileAttributes;
import java.util.HashMap;
import java.util.Map;

/**
 * Utilities to help transforming BasicFileAttributes to Map
 */
public abstract class AttributesUtils {

    /**
     * Given a BasicFileAttributes not null then return a Map
     * with the keys as the fields of the BasicFileAttributes or PosixFileAttributes and the values
     * with the content of the fields
     *
     * @param attr BasicFileAttributes
     * @return Map String Object never null
     */
    public static Map<String, Object> fileAttributeToMap(BasicFileAttributes attr) {
        Map<String, Object> result = new HashMap<>();
        result.put("creationTime", attr.creationTime());
        result.put("fileKey", attr.fileKey());
        result.put("isDirectory", attr.isDirectory());
        result.put("isOther", attr.isOther());
        result.put("isRegularFile", attr.isRegularFile());
        result.put("isSymbolicLink", attr.isSymbolicLink());
        result.put("lastAccessTime", attr.lastAccessTime());
        result.put("lastModifiedTime", attr.lastModifiedTime());
        result.put("size", attr.size());

        if (attr instanceof PosixFileAttributes) {
            PosixFileAttributes posixAttr = (PosixFileAttributes) attr;
            result.put("permissions", posixAttr.permissions());
            result.put("owner", posixAttr.owner());
            result.put("group", posixAttr.group());
        }

        return result;
    }

    /**
     * transform the java.nio.file.attribute.BasicFileAttributes to Map filtering by the keys
     * given in the filters param
     *
     * @param attr    BasicFileAttributes not null to tranform to map
     * @param filters String[] filters
     * @return Map String Object with the same keys as the filters
     */
    public static Map<String, Object> fileAttributeToMap(BasicFileAttributes attr, String[] filters) {
        Map<String, Object> result = new HashMap<>();

        for (String filter : filters) {
            filter = filter.replace("basic:", "");
            filter = filter.replace("posix:", "");
            switch (filter) {
                case "creationTime":
                    result.put("creationTime", attr.creationTime());
                    break;
                case "fileKey":
                    result.put("fileKey", attr.fileKey());
                    break;
                case "isDirectory":
                    result.put("isDirectory", attr.isDirectory());
                    break;
                case "isOther":
                    result.put("isOther", attr.isOther());
                    break;
                case "isRegularFile":
                    result.put("isRegularFile", attr.isRegularFile());
                    break;
                case "isSymbolicLink":
                    result.put("isSymbolicLink", attr.isSymbolicLink());
                    break;
                case "lastAccessTime":
                    result.put("lastAccessTime", attr.lastAccessTime());
                    break;
                case "lastModifiedTime":
                    result.put("lastModifiedTime", attr.lastModifiedTime());
                    break;
                case "size":
                    result.put("size", attr.size());
                    break;
                case "permissions":
                    result.put("permissions", ((PosixFileAttributes)attr).permissions());
                    break;
                case "group":
                    result.put("group", ((PosixFileAttributes)attr).group());
                    break;
                case "owner":
                    result.put("owner", ((PosixFileAttributes)attr).owner());
                    break;
            }
        }

        return result;
    }
}
