package cromwell.engine.io.nio

import java.io._
import java.nio.charset.StandardCharsets

import akka.stream.scaladsl.Flow
import cats.effect.{IO, Timer}
import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import common.util.IORetry
import cromwell.core.io._
import cromwell.core.path.Path
import cromwell.engine.io.IoActor._
import cromwell.engine.io.RetryableRequestSupport.{isRetryable, isInfinitelyRetryable}
import cromwell.engine.io.{IoAttempts, IoCommandContext}
import cromwell.filesystems.drs.DrsPath
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.oss.OssPath
import cromwell.filesystems.s3.S3Path
import cromwell.util.TryWithResource._

import scala.concurrent.ExecutionContext
object NioFlow {
  val NoopOnRetry: IoCommandContext[_] => Throwable => Unit = _ => _ => ()
}

/**
  * Flow that executes IO operations by calling java.nio.Path methods
  */
class NioFlow(parallelism: Int,
              onRetryCallback: IoCommandContext[_] => Throwable => Unit = NioFlow.NoopOnRetry,
              nbAttempts: Int = MaxAttemptsNumber)(implicit ec: ExecutionContext) {
  
  implicit private val timer: Timer[IO] = IO.timer(ec)
  
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
      isRetryable = isRetryable,
      isInfinitelyRetryable = isInfinitelyRetryable,
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

  private[io] val flow =
    Flow[DefaultCommandContext[_]].mapAsyncUnordered[IoResult](parallelism)(processCommand.andThen(_.unsafeToFuture()))

  private def copy(copy: IoCopyCommand) = IO {
    createDirectories(copy.destination)
    copy.source.copyTo(copy.destination, overwrite = true)
    ()
  }

  private def write(write: IoWriteCommand) = IO {
    createDirectories(write.file)
    write.file.writeContent(write.content)(write.openOptions, StandardCharsets.UTF_8, write.compressPayload)
    ()
  }

  private def delete(delete: IoDeleteCommand) = IO {
    delete.file.delete(delete.swallowIOExceptions)
    ()
  }

  private def readAsString(read: IoContentAsStringCommand) = IO {
    new String(
      read.file.limitFileContent(read.options.maxBytes, read.options.failOnOverflow),
      StandardCharsets.UTF_8
    )
      .replaceAll("\\r\\n", "\\\n")
  }

  private def size(size: IoSizeCommand) = IO {
    size.file.size
  }

  private def hash(hash: IoHashCommand): IO[String] = {
    hash.file match {
      case gcsPath: GcsPath => IO.fromTry { gcsPath.objectBlobId.map(gcsPath.cloudStorage.get(_).getCrc32c) }
      case drsPath: DrsPath => getFileHashForDrsPath(drsPath)
      case s3Path: S3Path => IO { s3Path.eTag }
      case ossPath: OssPath => IO { ossPath.eTag}
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
    exists.file.withReader { reader =>
      Stream.continually(reader.readLine()).takeWhile(_ != null).toList
    }
  }

  private def isDirectory(isDirectory: IoIsDirectoryCommand) = IO {
    isDirectory.file.isDirectory
  }

  private def createDirectories(path: Path) = path.parent.createDirectories()

  private def getFileHashForDrsPath(drsPath: DrsPath): IO[String] = {
    val drsFileSystemProvider = drsPath.drsPath.getFileSystem.provider.asInstanceOf[DrsCloudNioFileSystemProvider]

    //Since this does not actually do any IO work it is not wrapped in IO
    val fileAttributesOption = drsFileSystemProvider.fileProvider.fileAttributes(drsPath.drsPath.cloudHost, drsPath.drsPath.cloudPath)

    fileAttributesOption match {
      case Some(fileAttributes) =>
        val fileHashIO = IO { fileAttributes.fileHash }

        fileHashIO.flatMap({
          case Some(fileHash) => IO.pure(fileHash)
          case None => IO.raiseError(new IOException(s"Error while resolving DRS path $drsPath. The response from Martha doesn't contain the 'md5' hash for the file."))
        })
      case None => IO.raiseError(new IOException(s"Error getting file hash of DRS path $drsPath. Reason: File attributes class DrsCloudNioRegularFileAttributes wasn't defined in DrsCloudNioFileProvider."))
    }
  }
}
