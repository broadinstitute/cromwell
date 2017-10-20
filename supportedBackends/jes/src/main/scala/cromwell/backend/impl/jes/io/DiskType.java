package cromwell.backend.impl.jes.io;

public enum DiskType {
    LOCAL("LOCAL", "LOCAL_SSD"),
    SSD("SSD", "PERSISTENT_SSD"),
    HDD("HDD", "PERSISTENT_HDD");

    public final String diskTypeName;
    public final String googleTypeName;

    DiskType(final String diskTypeName, final String googleTypeName) {
        this.diskTypeName = diskTypeName;
        this.googleTypeName = googleTypeName;
    }
}
