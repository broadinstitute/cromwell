package cloud.nio.spi

import java.nio.file.attribute.{BasicFileAttributes, FileTime}

import cloud.nio.spi.CloudNioFileAttributes._

trait CloudNioFileAttributes extends BasicFileAttributes {
  def fileHash: Option[FileHash]
}

object CloudNioFileAttributes {
  val FileTimeZero: FileTime = FileTime.fromMillis(0)
}

case class FileHash(hashType: FileHash.Value, hash: String)

object FileHash extends Enumeration {
  val Crc32c: FileHash.Value = Value("crc32c")
  val Etag: FileHash.Value = Value("etag")
  val Md5: FileHash.Value = Value("md5")
  val Sha256: FileHash.Value = Value("sha256")
}

sealed trait ChecksumResult
case class ChecksumSuccess() extends ChecksumResult
case class ChecksumFailure(calculatedHash: String) extends ChecksumResult

trait CloudNioRegularFileAttributes extends CloudNioFileAttributes {
  override def lastAccessTime(): FileTime = FileTimeZero

  override def creationTime(): FileTime = lastModifiedTime()

  final override val isRegularFile: Boolean = true

  final override val isDirectory: Boolean = false

  final override val isSymbolicLink: Boolean = false

  final override val isOther: Boolean = false
}

final case class CloudNioDirectoryAttributes(path: CloudNioPath) extends CloudNioFileAttributes {
  override val lastModifiedTime: FileTime = FileTimeZero

  override val lastAccessTime: FileTime = FileTimeZero

  override val creationTime: FileTime = FileTimeZero

  override val isRegularFile: Boolean = false

  override val isDirectory: Boolean = true

  override val isSymbolicLink: Boolean = false

  override val isOther: Boolean = false

  override val size: Long = 0

  override val fileKey: AnyRef = path

  override val fileHash: Option[FileHash] = None
}
