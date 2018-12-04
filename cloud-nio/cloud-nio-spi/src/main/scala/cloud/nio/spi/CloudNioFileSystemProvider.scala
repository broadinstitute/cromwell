package cloud.nio.spi

import java.net.URI
import java.nio.channels.SeekableByteChannel
import java.nio.file._
import java.nio.file.attribute.{BasicFileAttributeView, BasicFileAttributes, FileAttribute, FileAttributeView}
import java.nio.file.spi.FileSystemProvider

import com.typesafe.config.{Config, ConfigFactory}

import scala.collection.JavaConverters._
import net.ceedubs.ficus.Ficus._

/**
  * Copy-port of https://github.com/broadinstitute/cromwell/blob/fa6eba28e2cf6020b24a56bcb7237acd9c74a4ac/filesystems/oss/src/main/scala/cromwell/filesystems/oss/nio/CloudNioFileSystemProvider.scala
  * that seems to be a copy-port of https://github.com/GoogleCloudPlatform/google-cloud-java/blob/fad70bdfdbc88e5c2ddf13be7085a7e9963f66c8/google-cloud-clients/google-cloud-contrib/google-cloud-nio/src/main/java/com/google/cloud/storage/contrib/nio/CloudStorageFileSystemProvider.java
  */
abstract class CloudNioFileSystemProvider extends FileSystemProvider {

  def config: Config
  
  def usePseudoDirectories: Boolean = true

  def fileProvider: CloudNioFileProvider

  /**
    * Returns true for exceptions that should never be retried.
    */
  def isFatal(exception: Exception): Boolean

  /**
    * Returns true for exceptions that should not count against the retry count.
    */
  def isTransient(exception: Exception): Boolean

  lazy val retry: CloudNioRetry = new CloudNioRetry(config) {
    override def isTransient(exception: Exception): Boolean = CloudNioFileSystemProvider.this.isTransient(exception)

    override def isFatal(exception: Exception): Boolean = CloudNioFileSystemProvider.this.isFatal(exception)
  }

  def newCloudNioFileSystem(uriAsString: String, config: Config): CloudNioFileSystem = {
    val host = getHost(uriAsString)
    new CloudNioFileSystem(this, host)
  }

  def getHost(uriAsString: String): String = {
    require(uriAsString.startsWith(getScheme + "://"), s"Scheme does not match $getScheme")

    //In some cases for a URI, the host name is null. For example, for DRS urls like 'dos://dg.123/123-123-123',
    //even though 'dg.123' is a valid host, somehow since it does not conform to URI's standards, uri.getHost returns null. In such
    //cases, authority is used instead of host.
    val uri = new URI(uriAsString)
    val host = uri.getHost
    val hostOrAuthority = if (host == null) uri.getAuthority else host
    require(!hostOrAuthority.isEmpty, s"Bucket/Host is empty")

    hostOrAuthority
  }

  override def newFileSystem(uri: URI, env: java.util.Map[String, _]): CloudNioFileSystem = {
    val config = ConfigFactory.parseMap(env)
    newCloudNioFileSystem(uri.toString, config)
  }

  override def getFileSystem(uri: URI): CloudNioFileSystem = {
    val config = CloudNioFileSystemProvider.defaultConfig(uri.getScheme)
    newCloudNioFileSystem(uri.toString, config)
  }

  override def getPath(uri: URI): CloudNioPath = {
    getFileSystem(uri).getPath(uri.getPath)
  }

  override def newByteChannel(
    path: Path,
    options: java.util.Set[_ <: OpenOption],
    attrs: FileAttribute[_]*
  ): SeekableByteChannel = {
    val cloudNioPath = CloudNioPath.checkPath(path)

    for (opt <- options.asScala) {
      opt match {
        case StandardOpenOption.READ | StandardOpenOption.WRITE | StandardOpenOption.SPARSE |
            StandardOpenOption.TRUNCATE_EXISTING | StandardOpenOption.CREATE | StandardOpenOption.CREATE_NEW =>
        /* ok */
        case StandardOpenOption.APPEND | StandardOpenOption.DELETE_ON_CLOSE | StandardOpenOption.DSYNC |
            StandardOpenOption.SYNC =>
          throw new UnsupportedOperationException(opt.toString)
      }
    }

    if (options.contains(StandardOpenOption.READ) && options.contains(StandardOpenOption.WRITE)) {
      throw new UnsupportedOperationException("Cannot open a READ+WRITE channel")
    } else if (options.contains(StandardOpenOption.WRITE)) {
      cloudNioWriteChannel(retry, cloudNioPath)
    } else {
      cloudNioReadChannel(retry, cloudNioPath)
    }
  }
  
  protected def cloudNioReadChannel(retry: CloudNioRetry, cloudNioPath: CloudNioPath): CloudNioReadChannel = new CloudNioReadChannel(fileProvider, retry, cloudNioPath)
  protected def cloudNioWriteChannel(retry: CloudNioRetry, cloudNioPath: CloudNioPath): CloudNioWriteChannel = new CloudNioWriteChannel(fileProvider, retry, cloudNioPath)

  override def createDirectory(dir: Path, attrs: FileAttribute[_]*): Unit = retry.from(() => {
    val cloudNioPath = CloudNioPath.checkPath(dir)
    fileProvider.createDirectory(cloudNioPath.cloudHost, cloudNioPath.cloudPath)
  })

  override def deleteIfExists(path: Path): Boolean = {
    val cloudNioPath = CloudNioPath.checkPath(path)

    if (checkDirectoryExists(cloudNioPath)) {
      val hasObjects = retry.from(
        () => fileProvider.existsPaths(cloudNioPath.cloudHost, cloudNioPath.cloudPath)
      )
      if (hasObjects) {
        throw new UnsupportedOperationException("Can not delete a non-empty directory")
      } else {
        true
      }
    } else {
      retry.from(
        () => fileProvider.deleteIfExists(cloudNioPath.cloudHost, cloudNioPath.cloudPath)
      )
    }
  }

  override def delete(path: Path): Unit = {
    if (!deleteIfExists(path)) {
      val cloudNioPath = CloudNioPath.checkPath(path)
      throw new NoSuchFileException(cloudNioPath.uriAsString)
    }
  }

  override def copy(source: Path, target: Path, options: CopyOption*): Unit = {
    val sourceCloudNioPath = CloudNioPath.checkPath(source)
    val targetCloudNioPath = CloudNioPath.checkPath(target)

    if (sourceCloudNioPath != targetCloudNioPath) {
      retry.from(
        () =>
          fileProvider.copy(
            sourceCloudNioPath.cloudHost,
            sourceCloudNioPath.cloudPath,
            targetCloudNioPath.cloudHost,
            targetCloudNioPath.cloudPath
          )
      )
    }
  }

  override def move(source: Path, target: Path, options: CopyOption*): Unit = {
    for (option <- options) {
      if (option == StandardCopyOption.ATOMIC_MOVE)
        throw new AtomicMoveNotSupportedException(null, null, "Atomic move unsupported")
    }
    copy(source, target, options: _*)
    delete(source)
    ()
  }

  override def isSameFile(path: Path, path2: Path): Boolean = {
    CloudNioPath.checkPath(path).equals(CloudNioPath.checkPath(path2))
  }

  override def isHidden(path: Path): Boolean = {
    CloudNioPath.checkPath(path)
    false
  }

  override def getFileStore(path: Path): FileStore = throw new UnsupportedOperationException

  override def checkAccess(path: Path, modes: AccessMode*): Unit = {
    if (modes.nonEmpty) {
      throw new UnsupportedOperationException("Checking access modes is not supported, only file existence.")
    }

    val cloudNioPath = CloudNioPath.checkPath(path)

    val exists = checkDirectoryExists(cloudNioPath) || retry.from(
      () => fileProvider.existsPath(cloudNioPath.cloudHost, cloudNioPath.cloudPath)
    )

    if (!exists) {
      throw new NoSuchFileException(cloudNioPath.uriAsString)
    }
  }

  def checkDirectoryExists(cloudNioPath: CloudNioPath): Boolean = {
    // Anything that "seems" like a directory exists. Otherwise see if the path with a "/" contains files on the cloud.
    (usePseudoDirectories && cloudNioPath.seemsLikeDirectory) || retry.from(
      () => fileProvider.existsPaths(cloudNioPath.cloudHost, cloudNioPath.cloudPath + "/")
    )
  }

  override def getFileAttributeView[V <: FileAttributeView](
    path: Path,
    viewType: Class[V],
    options: LinkOption*
  ): V = {
    if (viewType != classOf[CloudNioFileAttributeView] && viewType != classOf[BasicFileAttributeView]) {
      throw new UnsupportedOperationException(viewType.getSimpleName)
    }

    val cloudNioPath = CloudNioPath.checkPath(path)
    val isDirectory = checkDirectoryExists(cloudNioPath)

    CloudNioFileAttributeView(fileProvider, retry, cloudNioPath, isDirectory).asInstanceOf[V]
  }

  override def readAttributes(path: Path, attributes: String, options: LinkOption*): java.util.Map[String, AnyRef] = {
    throw new UnsupportedOperationException
  }

  override def readAttributes[A <: BasicFileAttributes](
    path: Path,
    attributesType: Class[A],
    options: LinkOption*
  ): A = {
    if (attributesType != classOf[CloudNioFileAttributes] && attributesType != classOf[BasicFileAttributes]) {
      throw new UnsupportedOperationException(attributesType.getSimpleName)
    }

    val cloudNioPath = CloudNioPath.checkPath(path)

    if (checkDirectoryExists(cloudNioPath)) {
      CloudNioDirectoryAttributes(cloudNioPath).asInstanceOf[A]
    } else {
      retry
        .from(
          () => fileProvider.fileAttributes(cloudNioPath.cloudHost, cloudNioPath.cloudPath)
        )
        .map(_.asInstanceOf[A])
        .getOrElse(throw new NoSuchFileException(cloudNioPath.uriAsString))
    }
  }

  override def newDirectoryStream(dir: Path, filter: DirectoryStream.Filter[_ >: Path]): DirectoryStream[Path] = {
    val cloudNioPath = CloudNioPath.checkPath(dir)
    new CloudNioDirectoryStream(fileProvider, retry, cloudNioPath, filter)
  }

  override def setAttribute(path: Path, attribute: String, value: scala.Any, options: LinkOption*): Unit = {
    throw new UnsupportedOperationException
  }

  def canEqual(other: Any): Boolean = other.isInstanceOf[CloudNioFileSystemProvider]

  override def equals(other: Any): Boolean = other match {
    case that: CloudNioFileSystemProvider =>
      (that canEqual this) &&
        config == that.config
    case _ => false
  }

  override def hashCode(): Int = {
    val state = List(config)
    state.map(_.hashCode()).foldLeft(0)((a, b) => 31 * a + b)
  }
}

object CloudNioFileSystemProvider {

  def defaultConfig(scheme: String): Config = {
    ConfigFactory.load.getOrElse(s"cloud.nio.default.$scheme", ConfigFactory.empty)
  }
}
