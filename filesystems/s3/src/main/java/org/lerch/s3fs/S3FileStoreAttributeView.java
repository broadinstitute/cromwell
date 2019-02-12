package org.lerch.s3fs;

import java.nio.file.attribute.FileStoreAttributeView;
import java.util.Date;

public class S3FileStoreAttributeView implements FileStoreAttributeView {

    public static final String ATTRIBUTE_VIEW_NAME = "S3FileStoreAttributeView";

    private Date creationDate;
    private String name;
    private String ownerId;
    private String ownerDisplayName;

    public static enum AttrID {
        creationDate, name, ownerId, ownerDisplayName
    }

    public S3FileStoreAttributeView(Date creationDate, String name, String ownerId, String ownerDisplayName) {
        this.creationDate = creationDate;
        this.name = name;
        this.ownerId = ownerId;
        this.ownerDisplayName = ownerDisplayName;
    }

    @Override
    public String name() {
        return ATTRIBUTE_VIEW_NAME;
    }

    public Object getAttribute(String attribute) {
        return getAttribute(AttrID.valueOf(attribute));
    }

    private Object getAttribute(AttrID attrID) {
        switch (attrID) {
            case creationDate:
                return creationDate;
            case ownerDisplayName:
                return ownerDisplayName;
            case ownerId:
                return ownerId;
            default:
                return name;
        }
    }
}