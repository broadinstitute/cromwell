package cromwell.engine.io.oss


import akka.actor.{ActorSystem, Scheduler}
import akka.stream.scaladsl.Flow
import cromwell.core.io.IoCommand
import cromwell.core.retry.Retry
import cromwell.engine.io.IoActor
import cromwell.engine.io.IoActor.{IoResult, OssCommandContext}
import cromwell.filesystems.oss._
import java.nio.file.{Files, Paths}

import better.files.File.OpenOptions

import scala.concurrent.{ExecutionContext, Future}

class OssFlow(parallelism: Int, scheduler: Scheduler, nbAttempts: Int = IoActor.MaxAttemptsNumber)(implicit ec: ExecutionContext, actorSystem: ActorSystem) {
  private val processCommand: OssCommandContext[_] => Future[IoResult] = commandContext => {
    val operationResult = Retry.withRetry(
      () => handleSingleCommand(commandContext.request),
      maxRetries = Option(nbAttempts),
      backoff = IoCommand.defaultBackoff,
      isTransient = IoActor.isTransient,
      isFatal = IoActor.isFatal
    )

    operationResult map { (_, commandContext) } recoverWith {
      case failure => Future.successful(commandContext.fail(failure))
    }
  }

  private [oss] def handleSingleCommand(ioSingleCommand: OssIoCommand[_]) = {
    ioSingleCommand match {
      case copyCommand: OssIoCopyCommand => copy(copyCommand) map copyCommand.success
      case writeCommand: OssIoWriteCommand => write(writeCommand) map writeCommand.success
      case deleteCommand: OssIoDeleteCommand => delete(deleteCommand) map deleteCommand.success
      case sizeCommand: OssIoSizeCommand => size(sizeCommand) map sizeCommand.success
      case readAsStringCommand: OssIoContentAsStringCommand => readAsString(readAsStringCommand) map readAsStringCommand.success
      case hashCommand: OssIoHashCommand => hash(hashCommand) map hashCommand.success
      case touchCommand: OssIoTouchCommand => touch(touchCommand) map touchCommand.success
      case _ => Future.failed(new NotImplementedError("Method not implemented"))
    }
  }

  val flow = Flow[OssCommandContext[_]].mapAsyncUnordered[IoResult](parallelism)(processCommand)

  private def copy(copy: OssIoCopyCommand) = Future {
    copy.source.copyTo(copy.destination, true)
  }

  private def write(write: OssIoWriteCommand) = Future {
    val file = write.file
    val content = write.content
    val options = write.openOptions

    options match {
      case Seq(UploadFileOption) =>
        val bytes = Files.readAllBytes(Paths.get(write.content))
        file.writeByteArray(bytes)(OpenOptions.default)
      case _ =>
        file.writeByteArray(content.getBytes)(OpenOptions.default)
    }
  }

  private def delete(delete: OssIoDeleteCommand) = Future {
    delete.file.delete(true)
  }

  private def size(size: OssIoSizeCommand) = Future {
    size.file.size
  }

  private def readAsString(command: OssIoContentAsStringCommand) = Future {
    command.file.contentAsString
  }

  private def hash(hash: OssIoHashCommand) = Future {
    hash.file.md5
  }

  private def touch(touch :OssIoTouchCommand) = Future {
    touch.file.touch()
  }

}
