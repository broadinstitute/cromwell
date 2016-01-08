package cromwell.engine.io.gcs

import java.io.{FileNotFoundException, OutputStream}
import java.net.URI
import java.nio.channels.{Channels, SeekableByteChannel}
import java.nio.file.DirectoryStream.Filter
import java.nio.file._
import java.nio.file.attribute.{BasicFileAttributes, FileAttribute, FileAttributeView}
import java.nio.file.spi.FileSystemProvider
import java.util
import java.util.Collections
import java.util.concurrent.{AbstractExecutorService, TimeUnit}

import com.google.api.services.storage.model.StorageObject
import com.google.cloud.hadoop.gcsio.{GoogleCloudStorageReadChannel, GoogleCloudStorageWriteChannel, ObjectWriteConditions}
import com.google.cloud.hadoop.util.{ApiErrorExtractor, AsyncWriteChannelOptions, ClientRequestHelper}

import scala.collection.JavaConverters._
import scala.concurrent.{ExecutionContext, ExecutionContextExecutorService}

object GcsFileSystemProvider {
  def instance(gcsInterface: GoogleCloudStorage)(implicit executionContext: ExecutionContext) = new GcsFileSystemProvider(gcsInterface, executionContext)
}

/**
  * Converts a Scala ExecutionContext to a Java ExecutorService.
  * https://groups.google.com/forum/#!topic/scala-user/ZyHrfzD7eX8
  */
object ExecutionContextExecutorServiceBridge {
  def apply(ec: ExecutionContext): ExecutionContextExecutorService = ec match {
    case null => throw new Throwable("Execution context cannot be null")
    case eces: ExecutionContextExecutorService => eces
    case executionContext => new AbstractExecutorService with ExecutionContextExecutorService {
      override def prepare(): ExecutionContext = executionContext
      override def isShutdown = false
      override def isTerminated = false
      override def shutdown() = ()
      override def shutdownNow() = Collections.emptyList[Runnable]
      override def execute(runnable: Runnable): Unit = executionContext execute runnable
      override def reportFailure(t: Throwable): Unit = executionContext reportFailure t
      override def awaitTermination(length: Long,unit: TimeUnit): Boolean = false
    }
  }
}

/**
  * Implements java.nio.FileSystemProvider for GoogleCloudStorage
  * This implementation is not complete and mostly a proof of concept that it's possible to *copy* around files from/to local/gcs.
  * Copying is the only functionality that has been successfully tested (same and cross filesystems).
  * @param googleCloudStorage must be properly set up (credentials) according to the context. Might be absorbed by this class eventually.
  * @param executionContext executionContext, will be used to perform async writes to GCS after being converted to a Java execution service
  */
class GcsFileSystemProvider private (googleCloudStorage: GoogleCloudStorage, executionContext: ExecutionContext) extends FileSystemProvider {

  private val executionService = ExecutionContextExecutorServiceBridge(executionContext)
  private val errorExtractor = new ApiErrorExtractor()

  private def checkExists(path: Path) = {
    if (!googleCloudStorage.exists(path.toString)) throw new FileNotFoundException(path.toString)
  }

  override def move(source: Path, target: Path, options: CopyOption*): Unit = throw new NotImplementedError()

  override def checkAccess(path: Path, modes: AccessMode*): Unit = {checkExists(path)}

  override def createDirectory(dir: Path, attrs: FileAttribute[_]*): Unit = {}

  override def getFileSystem(uri: URI): FileSystem = throw new UnsupportedOperationException()

  /**
    * Note: options and attributes are not honored.
    */
  override def newByteChannel(path: Path, options: util.Set[_ <: OpenOption], attrs: FileAttribute[_]*): SeekableByteChannel = {
    path match {
      case gcsPath: NioGcsPath =>
        new GoogleCloudStorageReadChannel(googleCloudStorage.client,
          gcsPath.bucket,
          gcsPath.objectName,
          errorExtractor,
          new ClientRequestHelper[StorageObject]()
        )
      case _ => throw new UnsupportedOperationException("Only Gcs paths are supported.")
    }
  }

  /**
    * Overrides the default implementation to provide a writable channel (which newByteChannel doesn't).
    * NOTE: options are not honored.
    */
  override def newOutputStream(path: Path, options: OpenOption*): OutputStream = {

    def initializeOutputStream(gcsPath: NioGcsPath) = {
      val channel = new GoogleCloudStorageWriteChannel(
        executionService,
        googleCloudStorage.client,
        new ClientRequestHelper[StorageObject](),
        gcsPath.bucket,
        gcsPath.objectName,
        AsyncWriteChannelOptions.newBuilder().build(),
        new ObjectWriteConditions(),
        Map.empty[String, String].asJava,
        "text/plain")
      channel.initialize()
      Channels.newOutputStream(channel)
    }

    path match {
      case gcsPath: NioGcsPath => initializeOutputStream(gcsPath)
      case _ => throw new UnsupportedOperationException("Only Gcs paths are supported.")
    }
  }

  override def isHidden(path: Path): Boolean = throw new NotImplementedError()

  override def copy(source: Path, target: Path, options: CopyOption*): Unit = {
    (source, target) match {
      case (s: NioGcsPath, d: NioGcsPath) => googleCloudStorage.copy(s, d)
      case _ => throw new UnsupportedOperationException(s"Can only copy from GCS to GCS: $source or $target is not a GCS path")
    }
  }

  override def delete(path: Path): Unit = {
    path match {
      case gcs: GcsPath => googleCloudStorage.deleteObject(gcs)
    }
  }

  override def newDirectoryStream(dir: Path, filter: Filter[_ >: Path]): DirectoryStream[Path] = throw new NotImplementedError()

  override def setAttribute(path: Path, attribute: String, value: scala.Any, options: LinkOption*): Unit = throw new NotImplementedError()

  override def getPath(uri: URI): Path = throw new NotImplementedError()

  override def newFileSystem(uri: URI, env: util.Map[String, _]): FileSystem = throw new NotImplementedError()

  override def readAttributes[A <: BasicFileAttributes](path: Path, `type`: Class[A], options: LinkOption*): A = {
    checkExists(path)
    new GcsFileAttributes(path).asInstanceOf[A]
  }

  override def readAttributes(path: Path, attributes: String, options: LinkOption*): util.Map[String, AnyRef] = throw new NotImplementedError()

  override def isSameFile(path: Path, path2: Path): Boolean = throw new NotImplementedError()

  override def getFileAttributeView[V <: FileAttributeView](path: Path, `type`: Class[V], options: LinkOption*): V = throw new NotImplementedError()

  override def getFileStore(path: Path): FileStore = throw new NotImplementedError()

  override def getScheme: String = GcsFileSystem.Protocol
}
