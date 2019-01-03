package cloud.nio.spi

import java.io.FileNotFoundException
import java.nio.file.attribute.{BasicFileAttributeView, FileTime}

final case class CloudNioFileAttributeView(
  fileProvider: CloudNioFileProvider,
  retry: CloudNioRetry,
  cloudNioPath: CloudNioPath,
  isDirectory: Boolean
) extends BasicFileAttributeView {
  override def name(): String = CloudNioFileAttributeView.Name

  override def readAttributes(): CloudNioFileAttributes = {
    if (isDirectory) {
      CloudNioDirectoryAttributes(cloudNioPath)
    } else {
      retry
        .from(
          () => fileProvider.fileAttributes(cloudNioPath.cloudHost, cloudNioPath.cloudPath)
        )
        .getOrElse(throw new FileNotFoundException(cloudNioPath.uriAsString))
    }
  }

  override def setTimes(lastModifiedTime: FileTime, lastAccessTime: FileTime, createTime: FileTime): Unit = {
    throw new UnsupportedOperationException
  }
}

object CloudNioFileAttributeView {
  val Name = "cloud"
}
