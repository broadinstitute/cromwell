package cromwell.core.path.proxy

import java.lang.Iterable
import java.nio.file._
import java.nio.file.attribute.UserPrincipalLookupService
import java.nio.file.spi.FileSystemProvider
import java.util

class FileSystemProxy(delegate: FileSystem, injectedProvider: FileSystemProvider) extends FileSystem {

  override def provider(): FileSystemProvider = injectedProvider

  /* delegated */
  override def supportedFileAttributeViews(): util.Set[String] = delegate.supportedFileAttributeViews()
  override def getSeparator: String = delegate.getSeparator
  override def getRootDirectories: Iterable[Path] = delegate.getRootDirectories
  override def newWatchService(): WatchService = delegate.newWatchService()
  override def getFileStores: Iterable[FileStore] = delegate.getFileStores
  override def isReadOnly: Boolean = delegate.isReadOnly
  override def getPath(first: String, more: String*): Path = new PathProxy(delegate.getPath(first, more: _*), this)
  override def isOpen: Boolean = delegate.isOpen
  override def close(): Unit = delegate.close()
  override def getPathMatcher(syntaxAndPattern: String): PathMatcher = delegate.getPathMatcher(syntaxAndPattern)
  override def getUserPrincipalLookupService: UserPrincipalLookupService = delegate.getUserPrincipalLookupService
}
