package cloud.nio.impl.ftp

import java.nio.file.attribute.FileTime

import cloud.nio.spi.CloudNioRegularFileAttributes
import org.apache.commons.net.ftp.FTPFile

class FtpCloudNioRegularFileAttributes(file: FTPFile, key: String) extends CloudNioRegularFileAttributes {
  override def fileHash = None
  override def lastModifiedTime() = FileTime.from(file.getTimestamp.toInstant)
  override def size() = file.getSize
  override def fileKey() = key
}
