package org.broadinstitute.manifestcreator.model;

public class ReferenceFile {
  private String path;
  private long crc32c;

  public String getPath() {
    return path;
  }

  public void setPath(String path) {
    this.path = path;
  }

  public long getCrc32c() {
    return crc32c;
  }

  public void setCrc32c(long crc32c) {
    this.crc32c = crc32c;
  }
}
