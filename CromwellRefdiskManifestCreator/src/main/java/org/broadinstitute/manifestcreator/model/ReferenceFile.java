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

  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (o == null || getClass() != o.getClass()) return false;

    ReferenceFile that = (ReferenceFile) o;

    if (getCrc32c() != that.getCrc32c()) return false;
    return getPath().equals(that.getPath());
  }

  @Override
  public int hashCode() {
    int result = getPath().hashCode();
    result = 31 * result + (int) (getCrc32c() ^ (getCrc32c() >>> 32));
    return result;
  }
}
