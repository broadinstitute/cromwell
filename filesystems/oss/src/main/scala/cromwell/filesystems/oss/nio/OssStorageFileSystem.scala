package cromwell.filesystems.oss.nio


import java.nio.file._
import java.nio.file.attribute.UserPrincipalLookupService
import java.util.Objects
import java.{lang, util}

import scala.collection.JavaConverters._


object OssStorageFileSystem {
  val SEPARATOR: String = "/"
  val URI_SCHEMA: String = "oss"
  val OSS_VIEW = "oss"
  val BASIC_VIEW = "basic"

  def apply(provider: OssStorageFileSystemProvider, bucket: String, config: OssStorageConfiguration): OssStorageFileSystem = {
    val res = new OssStorageFileSystem(bucket, config)
    res.internalProvider = provider
    res
  }
}



case class OssStorageFileSystem(bucket: String, config: OssStorageConfiguration) extends FileSystem {

  var internalProvider: OssStorageFileSystemProvider = OssStorageFileSystemProvider(config)

  override def provider: OssStorageFileSystemProvider = internalProvider

  override def getPath(first: String, more: String*): OssStoragePath = OssStoragePath.getPath(this, first, more: _*)

  override def close(): Unit = {
    // do nothing currently.
  }

  override def isOpen: Boolean = {
    true
  }

  override def isReadOnly: Boolean = {
    false
  }

  override def getSeparator: String = {
    OssStorageFileSystem.SEPARATOR
  }

  override def getRootDirectories: lang.Iterable[Path] = {
    Set[Path](OssStoragePath.getPath(this, UnixPath.ROOT_PATH)).asJava
  }

  override def getFileStores: lang.Iterable[FileStore] = {
    Set.empty[FileStore].asJava
  }

  override def getPathMatcher(syntaxAndPattern: String): PathMatcher = {
    FileSystems.getDefault.getPathMatcher(syntaxAndPattern)
  }

  override def getUserPrincipalLookupService: UserPrincipalLookupService = {
    throw new UnsupportedOperationException()
  }

  override def newWatchService(): WatchService = {
    throw new UnsupportedOperationException()
  }

  override def supportedFileAttributeViews(): util.Set[String] = {
    Set(OssStorageFileSystem.OSS_VIEW, OssStorageFileSystem.BASIC_VIEW).asJava
  }

  override def equals(obj: scala.Any): Boolean = {
    this == obj ||
      obj.isInstanceOf[OssStorageFileSystem] &&
        obj.asInstanceOf[OssStorageFileSystem].config.equals(config) &&
        obj.asInstanceOf[OssStorageFileSystem].bucket.equals(bucket)
  }

  override def hashCode(): Int = Objects.hash(bucket)

  override def toString: String = OssStorageFileSystem.URI_SCHEMA + "://" + bucket
}

