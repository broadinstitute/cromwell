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
  def getInstance(gcsInterface: GoogleCloudStorage)(implicit executionContext: ExecutionContext) = new GcsFileSystemProvider(gcsInterface, executionContext)
}

/**
  * Converts a Scala ExecutionContext to a Java ExecutorService.
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

class GcsFileSystemProvider(googleCloudStorage: GoogleCloudStorage, executionContext: ExecutionContext) extends FileSystemProvider {

  private val executionService = ExecutionContextExecutorServiceBridge(executionContext)
  private val errorExtractor = new ApiErrorExtractor()

  private def checkExists(path: Path) = {
    if(!googleCloudStorage.exists(path.toString)) throw new FileNotFoundException(path.toString)
  }

  override def move(source: Path, target: Path, options: CopyOption*): Unit = ???

  override def checkAccess(path: Path, modes: AccessMode*): Unit = {checkExists(path)}

  override def createDirectory(dir: Path, attrs: FileAttribute[_]*): Unit = {}

  override def getFileSystem(uri: URI): FileSystem = ???

  override def newByteChannel(path: Path, options: util.Set[_ <: OpenOption], attrs: FileAttribute[_]*): SeekableByteChannel = {
    path match {
      case gcsPath: GcsPath =>
        new GoogleCloudStorageReadChannel(googleCloudStorage.client,
          gcsPath.bucket,
          gcsPath.objectName,
          errorExtractor,
          new ClientRequestHelper[StorageObject]()
        )
      case _ => throw new UnsupportedOperationException("Only Gcs paths are supported.")
    }
  }

  override def newOutputStream(path: Path, options: OpenOption*): OutputStream = {
    path match {
      case gcsPath: GcsPath =>
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
      case _ => throw new UnsupportedOperationException("Only Gcs paths are supported.")
    }
  }

  override def isHidden(path: Path): Boolean = ???

  override def copy(source: Path, target: Path, options: CopyOption*): Unit = {
    (source, target) match {
      case (s: GcsPath, d: GcsPath) => googleCloudStorage.copy(s, d)
      case _ => throw new UnsupportedOperationException("Can only copy from GCS to GCS")
    }
  }

  override def delete(path: Path): Unit = {
    path match {
      case gcs: GcsPath => googleCloudStorage.deleteObject(gcs)
    }
  }

  override def newDirectoryStream(dir: Path, filter: Filter[_ >: Path]): DirectoryStream[Path] = ???

  override def setAttribute(path: Path, attribute: String, value: scala.Any, options: LinkOption*): Unit = ???

  override def getPath(uri: URI): Path = ???

  override def newFileSystem(uri: URI, env: util.Map[String, _]): FileSystem = ???

  override def readAttributes[A <: BasicFileAttributes](path: Path, `type`: Class[A], options: LinkOption*): A = {
    checkExists(path)
    new GcsFileAttributes(path).asInstanceOf[A]
  }

  override def readAttributes(path: Path, attributes: String, options: LinkOption*): util.Map[String, AnyRef] = ???

  override def isSameFile(path: Path, path2: Path): Boolean = ???

  override def getFileAttributeView[V <: FileAttributeView](path: Path, `type`: Class[V], options: LinkOption*): V = ???

  override def getFileStore(path: Path): FileStore = ???

  override def getScheme: String = ???
}
