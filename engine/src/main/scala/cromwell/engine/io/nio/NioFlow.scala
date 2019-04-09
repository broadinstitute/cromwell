package cromwell.engine.io.nio

import java.io._
import java.nio.charset.StandardCharsets

import akka.actor.Scheduler
import akka.stream.scaladsl.Flow
import cats.effect.IO
import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import common.util.IORetry
import cromwell.core.io._
import cromwell.core.path.Path
import cromwell.engine.io.IoActor._
import cromwell.engine.io.{IoAttempts, IoCommandContext}
import cromwell.filesystems.drs.DrsPath
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.s3.S3Path
import cromwell.util.TryWithResource._

import scala.concurrent.ExecutionContext
import scala.io.Codec
import scala.util.Failure
object NioFlow {
  def NoopOnRetry(context: IoCommandContext[_])(failure: Throwable) = ()
}

/**
  * Flow that executes IO operations by calling java.nio.Path methods
  */
class NioFlow(parallelism: Int,
              scheduler: Scheduler,
              onRetryCallback: IoCommandContext[_] => Throwable => Unit = NioFlow.NoopOnRetry,
              nbAttempts: Int = MaxAttemptsNumber)(implicit ec: ExecutionContext) {
  
  implicit private val timer = IO.timer(ec)
  
  private val processCommand: DefaultCommandContext[_] => IO[IoResult] = commandContext => {

    val onRetry: (Throwable, IoAttempts) => IoAttempts = (t, s) => {
      onRetryCallback(commandContext)
      IoAttempts.updateState(t, s)
    }

    val operationResult = IORetry.withRetry(
      handleSingleCommand(commandContext.request),
      IoAttempts(1),
      maxRetries = Option(nbAttempts),
      backoff = IoCommand.defaultBackoff,
      isTransient = isTransient,
      isFatal = isFatal,
      onRetry = onRetry
    )

    val io = for {
      _ <- IO.shift(ec)
      result <- operationResult
    } yield (result, commandContext)
    
     io handleErrorWith {
      failure => IO.pure(commandContext.fail(failure))
    }
  }

  private [nio] def handleSingleCommand(ioSingleCommand: IoCommand[_]): IO[IoSuccess[_]] = {
    ioSingleCommand match {
      case copyCommand: IoCopyCommand => copy(copyCommand) map copyCommand.success
      case writeCommand: IoWriteCommand => write(writeCommand) map writeCommand.success
      case deleteCommand: IoDeleteCommand => delete(deleteCommand) map deleteCommand.success
      case sizeCommand: IoSizeCommand => size(sizeCommand) map sizeCommand.success
      case readAsStringCommand: IoContentAsStringCommand => readAsString(readAsStringCommand) map readAsStringCommand.success
      case hashCommand: IoHashCommand => hash(hashCommand) map hashCommand.success
      case touchCommand: IoTouchCommand => touch(touchCommand) map touchCommand.success
      case existsCommand: IoExistsCommand => exists(existsCommand) map existsCommand.success
      case readLinesCommand: IoReadLinesCommand => readLines(readLinesCommand) map readLinesCommand.success
      case isDirectoryCommand: IoIsDirectoryCommand => isDirectory(isDirectoryCommand) map isDirectoryCommand.success
      case _ => IO.raiseError(new UnsupportedOperationException("Method not implemented"))
    }
  }

  val flow = Flow[DefaultCommandContext[_]].mapAsyncUnordered[IoResult](parallelism)(processCommand.andThen(_.unsafeToFuture()))

  private def copy(copy: IoCopyCommand) = IO {
    createDirectories(copy.destination)
    copy.source.copyTo(copy.destination, copy.overwrite)
    ()
  }

  private def write(write: IoWriteCommand) = IO {
    createDirectories(write.file)
    write.file.writeContent(write.content)(write.openOptions, Codec.UTF8)
    ()
  }

  private def delete(delete: IoDeleteCommand) = IO {
    delete.file.delete(delete.swallowIOExceptions)
    ()
  }

  private def readAsString(read: IoContentAsStringCommand) = IO {
    new String(
      limitFileContent(read.file, read.options.maxBytes, read.options.failOnOverflow),
      StandardCharsets.UTF_8
    )
  }

  private def size(size: IoSizeCommand) = IO {
    size.file.size
  }

  private def hash(hash: IoHashCommand): IO[String] = {
    hash.file match {
      case gcsPath: GcsPath => IO { gcsPath.cloudStorage.get(gcsPath.blob).getCrc32c }
      case drsPath: DrsPath => getFileHashForDrsPath(drsPath)
      case s3Path: S3Path => IO { s3Path.eTag }
      case path =>
        IO.fromEither(
        tryWithResource(() => path.newInputStream) { inputStream =>
          org.apache.commons.codec.digest.DigestUtils.md5Hex(inputStream)
        }.toEither
      )
    }
  }

  private def touch(touch: IoTouchCommand) = IO {
    touch.file.touch()
  }

  private def exists(exists: IoExistsCommand) = IO {
    exists.file.exists
  }

  private def readLines(exists: IoReadLinesCommand) = IO {
    withReader(exists.file) { reader =>
      Stream.continually(reader.readLine()).takeWhile(_ != null).toList
    }
  }

  private def isDirectory(isDirectory: IoIsDirectoryCommand) = IO {
    isDirectory.file.isDirectory
  }

  private def createDirectories(path: Path) = path.parent.createDirectories()

  /*
   * The input stream will be closed when this method returns, which means the f function
   * cannot leak an open stream.
   */
  private def withReader[A](file: Path)(f: BufferedReader => A): A = {
    
    // Use an input reader to convert the byte stream to character stream. Buffered reader for efficiency.
    tryWithResource(() => new BufferedReader(new InputStreamReader(file.mediaInputStream, Codec.UTF8.name)))(f).recoverWith({
      case failure => Failure(new IOException(s"Could not read from ${file.pathAsString}: ${failure.getMessage}", failure))
    }).get
  }

  /**
    * Returns an Array[Byte] from a Path. Limit the array size to "limit" byte if defined.
    * @throws IOException if failOnOverflow is true and the file is larger than limit
    */
  private def limitFileContent(file: Path, limit: Option[Int], failOnOverflow: Boolean) = withReader(file) { reader =>
    val bytesIterator = Iterator.continually(reader.read).takeWhile(_ != -1).map(_.toByte)
    // Take 1 more than the limit so that we can look at the size and know if it's overflowing
    val bytesArray = limit.map(l => bytesIterator.take(l + 1)).getOrElse(bytesIterator).toArray

    limit match {
      case Some(l) if failOnOverflow && bytesArray.length > l =>
        throw new IOException(s"File $file is larger than $l Bytes. Maximum read limits can be adjusted in the configuration under system.input-read-limits.")
      case Some(l) => bytesArray.take(l)
      case _ => bytesArray
    }
  }

  private def getFileHashForDrsPath(drsPath: DrsPath): IO[String] = {
    val drsFileSystemProvider = drsPath.drsPath.getFileSystem.provider.asInstanceOf[DrsCloudNioFileSystemProvider]

    //Since this does not actually do any IO work it is not wrapped in IO
    val fileAttributesOption = drsFileSystemProvider.fileProvider.fileAttributes(drsPath.drsPath.cloudHost, drsPath.drsPath.cloudPath)

    fileAttributesOption match {
      case Some(fileAttributes) => {
        val fileHashIO = IO { fileAttributes.fileHash }

        fileHashIO.flatMap({
          case Some(fileHash) => IO.pure(fileHash)
          case None => IO.raiseError(new IOException(s"Error while resolving DRS path $drsPath. The response from Martha doesn't contain the 'md5' hash for the file."))
        })
      }
      case None => IO.raiseError(new IOException(s"Error getting file hash of DRS path $drsPath. Reason: File attributes class DrsCloudNioRegularFileAttributes wasn't defined in DrsCloudNioFileProvider."))
    }
  }
}
