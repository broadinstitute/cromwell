package cromwell.engine.io.nio

import akka.actor.{ActorSystem, Scheduler}
import akka.stream.scaladsl.Flow
import cromwell.core.io._
import cromwell.core.retry.Retry
import cromwell.engine.io.IoActor
import cromwell.engine.io.IoActor.{DefaultCommandContext, IoResult}
import cromwell.filesystems.gcs.GcsPath
import cromwell.util.TryWithResource._

import scala.concurrent.{ExecutionContext, Future}
import scala.io.Codec

/**
  * Flow that executes IO operations by calling java.nio.Path methods
  */
class NioFlow(parallelism: Int, scheduler: Scheduler)(implicit ec: ExecutionContext, actorSystem: ActorSystem) {
  private val processCommand: DefaultCommandContext[_] => Future[IoResult] = commandContext => {
    val operationResult = Retry.withRetry(
      () => handleSingleCommand(commandContext.request),
      maxRetries = Option(3),
      backoff = IoCommand.defaultBackoff,
      isTransient = IoActor.isTransient,
      isFatal = IoActor.isFatal
    )
    
    operationResult map { (_, commandContext) } recoverWith {
      case failure => Future.successful(commandContext.fail(failure))
    }
  }
  
  private def handleSingleCommand(ioSingleCommand: IoCommand[_]) = {
    ioSingleCommand match {
      case copyCommand: IoCopyCommand => copy(copyCommand) map copyCommand.success
      case writeCommand: IoWriteCommand => write(writeCommand) map writeCommand.success
      case deleteCommand: IoDeleteCommand => delete(deleteCommand) map deleteCommand.success
      case sizeCommand: IoSizeCommand => size(sizeCommand) map sizeCommand.success
      case readAsStringCommand: IoContentAsStringCommand => readAsString(readAsStringCommand) map readAsStringCommand.success
      case hashCommand: IoHashCommand => hash(hashCommand) map hashCommand.success
      case _ => Future.failed(new NotImplementedError("Method not implemented"))
    }
  }
  
  val flow = Flow[DefaultCommandContext[_]].mapAsyncUnordered[IoResult](parallelism)(processCommand)
  
  private def copy(copy: IoCopyCommand) = Future {
    copy.source.copyTo(copy.destination, copy.overwrite)
    ()
  }

  private def write(write: IoWriteCommand) = Future {
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
}
