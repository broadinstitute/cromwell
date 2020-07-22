package org.broadinstitute.manifestcreator.model;

public class ReferenceFile {
  private String path;
  private long crc32;

  public String getPath() {
    return path;
  }

  public void setPath(String path) {
    this.path = path;
  }

  public long getCrc32() {
    return crc32;
  }

  public void setCrc32(long crc32) {
    this.crc32 = crc32;
  }
}
