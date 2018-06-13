package cromwell.backend.google.pipelines.common.io;

public enum DiskType {
    LOCAL("LOCAL", "local-ssd"),
    SSD("SSD", "pd-ssd"),
    HDD("HDD", "pd-standard");

    public final String diskTypeName;
    public final String googleTypeName;

    DiskType(final String diskTypeName, final String googleTypeName) {
        this.diskTypeName = diskTypeName;
        this.googleTypeName = googleTypeName;
    }
}
