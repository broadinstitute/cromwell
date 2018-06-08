package cromwell.engine.io.nio

import java.io.{IOException, InputStream}
import java.nio.charset.StandardCharsets

import akka.actor.{ActorSystem, Scheduler}
import akka.stream.scaladsl.Flow
import cromwell.core.io._
import cromwell.core.path.{DefaultPath, Path}
import cromwell.core.retry.Retry
import cromwell.engine.io.IoActor._
import cromwell.engine.io.IoCommandContext
import cromwell.filesystems.gcs.GcsPath
import cromwell.util.TryWithResource._

import scala.concurrent.{ExecutionContext, Future}
import scala.io.Codec

object NioFlow {
  def NoopOnRetry(context: IoCommandContext[_])(failure: Throwable) = ()
}

/**
  * Flow that executes IO operations by calling java.nio.Path methods
  */
class NioFlow(parallelism: Int,
              scheduler: Scheduler,
              onRetry: IoCommandContext[_] => Throwable => Unit = NioFlow.NoopOnRetry,
              nbAttempts: Int = MaxAttemptsNumber)(implicit ec: ExecutionContext, actorSystem: ActorSystem) {
  private val processCommand: DefaultCommandContext[_] => Future[IoResult] = commandContext => {
    val operationResult = Retry.withRetry(
      () => handleSingleCommand(commandContext.request),
      maxRetries = Option(nbAttempts),
      backoff = IoCommand.defaultBackoff,
      isTransient = isTransient,
      isFatal = isFatal,
      onRetry = onRetry(commandContext)
    )

    operationResult map { (_, commandContext) } recoverWith {
      case failure => Future.successful(commandContext.fail(failure))
    }
  }
  
  private [nio] def handleSingleCommand(ioSingleCommand: IoCommand[_]) = {
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
      case _ => Future.failed(new NotImplementedError("Method not implemented"))
    }
  }

  val flow = Flow[DefaultCommandContext[_]].mapAsyncUnordered[IoResult](parallelism)(processCommand)

  private def copy(copy: IoCopyCommand) = Future {
    createDirectoriesForSFSPath(copy.destination)
    copy.source.copyTo(copy.destination, copy.overwrite)
    ()
  }

  private def write(write: IoWriteCommand) = Future {
    createDirectoriesForSFSPath(write.file)
    write.file.write(write.content)(write.openOptions, Codec.UTF8)
    ()
  }

  private def delete(delete: IoDeleteCommand) = Future {
    delete.file.delete(delete.swallowIOExceptions)
    ()
  }

  private def readBytes(read: IoContentAsStringCommand, limit: Int, inputStream: InputStream) : String = {
    // Take 1 more than the limit so that we can look at the size and know if it's overflowing
    val bytesArray = Iterator.continually(inputStream.read).takeWhile(_ != -1).map(_.toByte).take(limit + 1).toArray
    if (read.options.failOnOverflow && bytesArray.length > limit)
      throw new IOException(s"File ${read.file.pathAsString} is larger than $limit Bytes. Maximum read limits can be adjusted in the configuration under system.input-read-limits.")
    else
      new String(bytesArray.take(limit), StandardCharsets.UTF_8)
  }

  private def readAsString(read: IoContentAsStringCommand) = {
    read.options.maxBytes match {
      case Some(limit) => Future(tryWithResource(() => read.file.mediaInputStream)(inputStream => readBytes(read, limit, inputStream)).get)
      case _ => Future(read.file.readContentAsString)
    }
  }

  private def size(size: IoSizeCommand) = Future {
    size.file.size
  }

  private def hash(hash: IoHashCommand) = {
    hash.file match {
      case gcsPath: GcsPath => Future { gcsPath.cloudStorage.get(gcsPath.blob).getCrc32c }
      case path => Future.fromTry(
        tryWithResource(() => path.newInputStream) { inputStream =>
          org.apache.commons.codec.digest.DigestUtils.md5Hex(inputStream)
        }
      )
    }
  }

  private def touch(touch: IoTouchCommand) = Future {
    touch.file.touch()
  }

  private def exists(exists: IoExistsCommand) = Future {
    exists.file.exists
  }

  private def readLines(exists: IoReadLinesCommand) = Future {
    exists.file.readAllLinesInFile
  }
  
  private def isDirectory(isDirectory: IoIsDirectoryCommand) = Future {
    isDirectory.file.isDirectory
  }

  private def createDirectoriesForSFSPath(path: Path) = path match {
    case _: DefaultPath => path.parent.createDirectories()
    case _ =>
  }
}
