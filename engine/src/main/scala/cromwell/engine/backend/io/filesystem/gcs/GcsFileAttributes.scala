package cromwell.engine.backend.io.filesystem.gcs

import java.nio.file.attribute.{BasicFileAttributes, FileTime}

import com.google.api.services.storage.Storage
import com.google.api.services.storage.model.StorageObject
import wdl4s.values._

class GcsFileAttributes(path: NioGcsPath, storageClient: Storage) extends BasicFileAttributes {
  override def fileKey(): AnyRef = path.toString.md5Sum
  override def isRegularFile: Boolean = throw new NotImplementedError("To be implemented when/if needed")
  override def isOther: Boolean = throw new NotImplementedError("To be implemented when/if needed")
  override def lastModifiedTime(): FileTime = throw new NotImplementedError("To be implemented when/if needed")
  override def size(): Long = {
    val getObject = storageClient.objects.get(path.bucket, path.objectName)
    val storageObject: StorageObject = getObject.execute()
    storageObject.getSize.longValue()
  }
  override def isDirectory: Boolean = path.isDirectory
  override def isSymbolicLink: Boolean = false
  override def creationTime(): FileTime = throw new NotImplementedError("To be implemented when/if needed")
  override def lastAccessTime(): FileTime = throw new NotImplementedError("To be implemented when/if needed")
}
