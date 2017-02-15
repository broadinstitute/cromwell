package cromwell.filesystems.gcs

import java.nio.file.{CopyOption, OpenOption}

import com.google.cloud.storage.Acl
import com.google.cloud.storage.contrib.nio.{CloudStorageOption, CloudStorageOptions}

sealed trait CromwellGcsOption extends OpenOption with CopyOption {
  def toGoogleOption: CloudStorageOption
}

case class CromwellGcsMimeType(mimeType: String) extends CromwellGcsOption {
  override def toGoogleOption: CloudStorageOption = CloudStorageOptions.withMimeType(mimeType)
}
case class CromwellGcsCacheControl(cacheControl: String) extends CromwellGcsOption {
  override def toGoogleOption: CloudStorageOption = CloudStorageOptions.withCacheControl(cacheControl)
}
case class CromwellGcsContentDisposition(contentDisposition: String) extends CromwellGcsOption {
  override def toGoogleOption: CloudStorageOption = CloudStorageOptions.withContentDisposition(contentDisposition)
}
case class CromwellGcsBlockSize(blockSize: Int) extends CromwellGcsOption {
  override def toGoogleOption: CloudStorageOption = CloudStorageOptions.withBlockSize(blockSize)
}
case class CromwellGcsContentEncoding(contentEncoding: String) extends CromwellGcsOption {
  override def toGoogleOption: CloudStorageOption = CloudStorageOptions.withContentEncoding(contentEncoding)
}
case class CromwellGcsUserMetadata(key: String, value: String) extends CromwellGcsOption {
  override def toGoogleOption: CloudStorageOption = CloudStorageOptions.withUserMetadata(key, value)
}
case class CromwellGcsAcl(acl: Acl) extends CromwellGcsOption {
  override def toGoogleOption: CloudStorageOption = CloudStorageOptions.withAcl(acl)
}
