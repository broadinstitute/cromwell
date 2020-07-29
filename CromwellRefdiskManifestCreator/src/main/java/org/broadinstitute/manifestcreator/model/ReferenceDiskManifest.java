package org.broadinstitute.manifestcreator.model;

import java.util.ArrayList;
import java.util.List;

public class ReferenceDiskManifest {
  private String imageIdentifier;
  private List<ReferenceFile> files;

  public String getImageIdentifier() {
    return imageIdentifier;
  }

  public void setImageIdentifier(String imageIdentifier) {
    this.imageIdentifier = imageIdentifier;
  }

  public List<ReferenceFile> getFiles() {
    if (files == null) {
      files = new ArrayList<>();
    }
    return files;
  }

  public void setFiles(List<ReferenceFile> files) {
    this.files = files;
  }
}
