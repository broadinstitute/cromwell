package cromwell.engine.io.nio

import akka.actor.{ActorSystem, Scheduler}
import akka.stream.scaladsl.Flow
import cromwell.core.io._
import cromwell.core.path.{DefaultPath, Path}
import cromwell.core.retry.Retry
import cromwell.engine.io.IoActor.{DefaultCommandContext, IoResult}
import cromwell.engine.io.{IoActor, IoCommandContext}
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
              nbAttempts: Int = IoActor.MaxAttemptsNumber)(implicit ec: ExecutionContext, actorSystem: ActorSystem) {
  private val processCommand: DefaultCommandContext[_] => Future[IoResult] = commandContext => {
    val operationResult = Retry.withRetry(
      () => handleSingleCommand(commandContext.request),
      maxRetries = Option(nbAttempts),
      backoff = IoCommand.defaultBackoff,
      isTransient = IoActor.isTransient,
      isFatal = IoActor.isFatal,
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

  private def readAsString(read: IoContentAsStringCommand) = Future {
    read.file.contentAsString
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

  private def createDirectoriesForSFSPath(path: Path) = path match {
    case _: DefaultPath => path.parent.createPermissionedDirectories()
    case _ =>
  }
}
