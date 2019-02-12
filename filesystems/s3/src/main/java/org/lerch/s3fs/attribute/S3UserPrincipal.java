package org.lerch.s3fs.attribute;

import java.nio.file.attribute.UserPrincipal;

public class S3UserPrincipal implements UserPrincipal {

    private String name;

    public S3UserPrincipal(String name) {
        this.name = name;
    }

    @Override
    public String getName() {
        return name;
    }
}
