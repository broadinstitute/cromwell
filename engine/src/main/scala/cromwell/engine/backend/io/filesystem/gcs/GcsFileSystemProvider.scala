package cromwell.engine.backend.io.filesystem.gcs

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

import com.google.api.client.googleapis.json.GoogleJsonResponseException
import com.google.api.services.storage.Storage
import com.google.api.services.storage.model.StorageObject
import com.google.cloud.hadoop.gcsio.{GoogleCloudStorageReadChannel, GoogleCloudStorageWriteChannel, ObjectWriteConditions}
import com.google.cloud.hadoop.util.{ApiErrorExtractor, AsyncWriteChannelOptions, ClientRequestHelper}

import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, ExecutionContextExecutorService}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object GcsFileSystemProvider {
  def apply()(implicit executionContext: ExecutionContext) = {
    new GcsFileSystemProvider(Failure(new Exception("No Storage object available")), executionContext)
  }
  def apply(storageClient: Storage)(implicit executionContext: ExecutionContext) = {
    new GcsFileSystemProvider(Success(storageClient), executionContext)
  }

  val defaultProvider = new GcsFileSystemProvider(Failure(new Exception("No Storage object available")), scala.concurrent.ExecutionContext.global)

  object AcceptAllFilter extends DirectoryStream.Filter[Path] {
    override def accept(entry: Path): Boolean = true
  }

  // To choose these numbers I first entered a prolonged period of personal consideration and deep thought.
  // Then, at the end of this time, I decided to just pick some numbers arbitrarily.
  private val retryInterval = 500 milliseconds
  private val retryCount = 3
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
  * @param storageClient Google API Storage object
  * @param executionContext executionContext, will be used to perform async writes to GCS after being converted to a Java execution service
  */
class GcsFileSystemProvider private (storageClient: Try[Storage], executionContext: ExecutionContext) extends FileSystemProvider {
  import GcsFileSystemProvider._

  // We want to throw an exception here if we try to use this class with a failed gcs interface
  lazy val client = storageClient.get
  private val executionService = ExecutionContextExecutorServiceBridge(executionContext)
  private val errorExtractor = new ApiErrorExtractor()
  def notAGcsPath(path: Path) = throw new IllegalArgumentException(s"$path is not a GCS path.")
  // Can't instantiate it now as it needs a GcsFileSystemProvider (this), which is not yet instantiated at this point..
  private var defaultFileSystem: Option[GcsFileSystem] = None

  private def exists(path: Path) = path match {
    case gcsPath: NioGcsPath =>
      val getObject = client.objects.get(gcsPath.bucket, gcsPath.objectName)

      def checkExists = {
        Try(getObject.execute) recoverWith {
          case ex: GoogleJsonResponseException if ex.getStatusCode == 404 =>
            if (gcsPath.isDirectory) Success(gcsPath) else Failure(ex)
        }
      }

      def tryExists(retries: Int): Boolean = checkExists match {
        case Success(_) => true
        case Failure(ex: GoogleJsonResponseException) if ex.getStatusCode == 404 && retries > 0 =>
          // FIXME remove this sleep
          Thread.sleep(retryInterval.toMillis.toInt)
          tryExists(retries - 1)
        case Failure(ex: GoogleJsonResponseException) if ex.getStatusCode == 404 => false
        case Failure(ex) => throw ex
      }

      if (!tryExists(retryCount)) throw new FileNotFoundException(path.toString)
    case _ => throw new FileNotFoundException(path.toString)
  }

  /**
    * Note: options and attributes are not honored.
    */
  override def newByteChannel(path: Path, options: util.Set[_ <: OpenOption], attrs: FileAttribute[_]*): SeekableByteChannel = {
    path match {
      case gcsPath: NioGcsPath =>
        new GoogleCloudStorageReadChannel(client,
          gcsPath.bucket,
          gcsPath.objectName,
          errorExtractor,
          new ClientRequestHelper[StorageObject]()
        )
      case _ => notAGcsPath(path)
    }
  }

  /**
    * Overrides the default implementation to provide a writable channel (which newByteChannel doesn't).
    * NOTE: options are not honored.
    */
  override def newOutputStream(path: Path, options: OpenOption*): OutputStream = {
    val contentType = options collectFirst {
      case e: ContentTypeOption.ContentType => e.toString
    } getOrElse ContentTypeOption.PlainText.toString

    def initializeOutputStream(gcsPath: NioGcsPath) = {
      val channel = new GoogleCloudStorageWriteChannel(
        executionService,
        client,
        new ClientRequestHelper[StorageObject](),
        gcsPath.bucket,
        gcsPath.objectName,
        AsyncWriteChannelOptions.newBuilder().build(),
        new ObjectWriteConditions(),
        Map.empty[String, String].asJava,
        contentType)
      channel.initialize()
      Channels.newOutputStream(channel)
    }

    path match {
      case gcsPath: NioGcsPath => initializeOutputStream(gcsPath)
      case _ => notAGcsPath(path)
    }
  }

  override def copy(source: Path, target: Path, options: CopyOption*): Unit = {
    (source, target) match {
      case (s: NioGcsPath, d: NioGcsPath) =>
        val storageObject = client.objects.get(s.bucket, s.objectName).execute
        client.objects.copy(s.bucket, s.objectName, d.bucket, d.objectName, storageObject).execute
      case _ => throw new UnsupportedOperationException(s"Can only copy from GCS to GCS: $source or $target is not a GCS path")
    }
  }

  override def delete(path: Path): Unit = {
    path match {
      case gcs: NioGcsPath => client.objects.delete(gcs.bucket, gcs.objectName).execute()
      case _ => notAGcsPath(path)
    }
  }

  override def readAttributes[A <: BasicFileAttributes](path: Path, `type`: Class[A], options: LinkOption*): A = path match {
    case gcsPath: NioGcsPath =>
      exists(path)
      new GcsFileAttributes(gcsPath, client).asInstanceOf[A]
    case _ => notAGcsPath(path)
  }

  override def move(source: Path, target: Path, options: CopyOption*): Unit = {
    (source, target) match {
      case (s: NioGcsPath, d: NioGcsPath) =>
        def moveInner = {
          val storageObject = client.objects.get(s.bucket, s.objectName).execute
          client.objects.rewrite(s.bucket, s.objectName, d.bucket, d.objectName, storageObject).execute
        }

        def tryMove(retries: Int): Unit = Try(moveInner) match {
          case Success(_) =>
          case Failure(ex: GoogleJsonResponseException) if ex.getStatusCode == 404 && retries > 0 =>
            /* Could not use TryUtil here because it requires a WorkflowLogger, which we can't get form here
             * TODO From a more general perspective we may need to add a logging capability to the IOInterface
             */
            Thread.sleep(retryInterval.toMillis)
            tryMove(retries - 1)
          case Failure(ex: GoogleJsonResponseException) if ex.getStatusCode == 404 =>
          case Failure(ex) => throw ex
        }

        tryMove(retryCount)

      case _ => throw new UnsupportedOperationException(s"Can only copy from GCS to GCS: $source or $target is not a GCS path")
    }
  }

  def crc32cHash(path: Path) = path match {
    case gcsDir: NioGcsPath =>
      val obj = client.objects().get(gcsDir.bucket, gcsDir.objectName).execute()
      obj.getCrc32c
    case _ => notAGcsPath(path)
  }

  override def checkAccess(path: Path, modes: AccessMode*): Unit = exists(path)
  override def createDirectory(dir: Path, attrs: FileAttribute[_]*): Unit = {}
  override def getFileSystem(uri: URI): FileSystem = getDefaultFileSystem

  def getDefaultFileSystem = defaultFileSystem match {
    case Some(fs) => fs
    case None =>
      val fs = GcsFileSystem(this)
      defaultFileSystem = Option(fs)
      fs
  }

  override def isHidden(path: Path): Boolean = throw new NotImplementedError()
  override def newDirectoryStream(dir: Path, filter: Filter[_ >: Path]): DirectoryStream[Path] = dir match {
    case gcsDir: NioGcsPath =>
      val listRequest = client.objects().list(gcsDir.bucket)
      listRequest.setPrefix(gcsDir.objectName)

      val paths = for {
        listedFile <- listRequest.execute().getItems.asScala
      } yield NioGcsPath(s"$getScheme${listedFile.getBucket}${GcsFileSystem.Separator}${listedFile.getName}")(dir.getFileSystem.asInstanceOf[GcsFileSystem])

      new DirectoryStream[Path] {
        override def iterator(): util.Iterator[Path] = paths.toIterator.asJava
        override def close(): Unit = {}
      }
    case _ => notAGcsPath(dir)
  }
  override def setAttribute(path: Path, attribute: String, value: scala.Any, options: LinkOption*): Unit = throw new NotImplementedError()
  override def getPath(uri: URI): Path = throw new NotImplementedError()
  override def newFileSystem(uri: URI, env: util.Map[String, _]): FileSystem = {
    throw new UnsupportedOperationException("GcsFileSystem provider doesn't support creation of new FileSystems at this time. Use getFileSystem instead.")
  }
  override def readAttributes(path: Path, attributes: String, options: LinkOption*): util.Map[String, AnyRef] = throw new NotImplementedError()
  override def isSameFile(path: Path, path2: Path): Boolean = throw new NotImplementedError()
  override def getFileAttributeView[V <: FileAttributeView](path: Path, `type`: Class[V], options: LinkOption*): V = throw new NotImplementedError()
  override def getFileStore(path: Path): FileStore = throw new NotImplementedError()
  override def getScheme: String = GcsFileSystem.Protocol
}
