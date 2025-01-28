package cloud.nio.impl.ftp

import java.nio.file.attribute.FileTime
import cloud.nio.spi.{CloudNioRegularFileAttributes, FileHash}
import org.apache.commons.net.ftp.FTPFile

class FtpCloudNioRegularFileAttributes(file: FTPFile, key: String) extends CloudNioRegularFileAttributes {
  override def fileHashes: List[FileHash] = List.empty
  override def lastModifiedTime() = FileTime.from(file.getTimestamp.toInstant)
  override def size() = file.getSize
  override def fileKey() = key
}
