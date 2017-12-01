package cromwell.engine.io.oss

import java.io.{BufferedReader, ByteArrayInputStream, InputStreamReader}

import scala.util.control._
import akka.actor.{ActorSystem, Scheduler}
import akka.stream.scaladsl.Flow
import cromwell.core.io.IoCommand
import cromwell.core.retry.Retry
import cromwell.engine.io.IoActor
import cromwell.engine.io.IoActor.{IoResult, OssCommandContext}
import cromwell.filesystems.oss._
import java.io.FileInputStream
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
    val src = copy.source
    val dest = copy.destination
    val ossClient = src.ossClient
    ossClient.copyObject(src.bucket, src.key, dest.bucket, dest.key)
  }

  private def write(write: OssIoWriteCommand) = Future {
    val file = write.file
    val content = write.content
    val ossClient = file.ossClient
    val options = write.openOptions

    options match {
      case Seq(UploadFileOption) =>
        val inputStream = new FileInputStream(content)
        ossClient.putObject(file.bucket, file.key, inputStream)
      case _ =>
        ossClient.putObject(file.bucket, file.key, new ByteArrayInputStream(content.getBytes()))
    }

  }

  private def delete(delete: OssIoDeleteCommand) = Future {
    val file = delete.file
    val ossClient = file.ossClient
    ossClient.deleteObject(file.bucket, file.key)
  }

  private def size(size: OssIoSizeCommand) = Future {
    val file = size.file
    val ossClient = file.ossClient
    val meta = ossClient.getObjectMetadata(file.bucket, file.key)
    meta.getContentLength
  }

  private def readAsString(command: OssIoContentAsStringCommand) = Future {
    val file = command.file
    val ossClient = file.ossClient
    val ossObject = ossClient.getObject(file.bucket, file.key)

    val reader = new BufferedReader(new InputStreamReader(ossObject.getObjectContent()))

    val ret = new StringBuilder
    val loop = new Breaks
    loop.breakable {
      while (true) {
        val line = reader.readLine()
        if (line == null) {
          loop.break
        }
        ret ++= (line + '\n')
      }
    }
    reader.close()
    ret.toString()
  }

  private def hash(hash: OssIoHashCommand) = Future {
    val file = hash.file
    val ossClient = file.ossClient
    val meta = ossClient.getObjectMetadata(file.bucket, file.key)
    meta.getContentMD5
  }

  private def touch(touch :OssIoTouchCommand) = Future {
    val file = touch.file
    val ossClient = file.ossClient
    val content: String = ""
    ossClient.putObject(file.bucket, file.key, new ByteArrayInputStream(content.getBytes()))
  }

}
