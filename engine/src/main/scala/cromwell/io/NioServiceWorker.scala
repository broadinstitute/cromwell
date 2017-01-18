package cromwell.io

import java.nio.file.Path

import akka.actor.{ActorLogging, Props}
import better.files.File
import better.files.File.OpenOptions
import cats.instances.try_._
import cats.syntax.functor._
import com.google.cloud.storage.StorageException
import cromwell.services.io._

import scala.concurrent.duration.FiniteDuration
import scala.util.{Failure, Success, Try}

class NioServiceWorker() extends IoServiceActor with ActorLogging {

  implicit val ec = context.dispatcher

  /* Command dispatch */
  override def handleCommand(command: IoActorCommand[_]) = {
    log.debug(s"processing IO command $command")
    
    val action = command match {
      case copyCommand : CopyCommand => 
        copy(copyCommand.source, copyCommand.destination, copyCommand.overwrite) map copyCommand.success
        
      case readAsStringCommand: ReadAsStringCommand => 
        readAsString(readAsStringCommand.file) map readAsStringCommand.success
        
      case sizeCommand: SizeCommand => 
        size(sizeCommand.file) map sizeCommand.success

      case writeCommand: WriteCommand => 
        write(writeCommand.file, writeCommand.content)(writeCommand.openOptions) map writeCommand.success

      case deleteCommand: DeleteCommand => 
        delete(deleteCommand.file, deleteCommand.swallowIOExceptions) map deleteCommand.success
    }
    
    action match {
      case Success(value) => value
      case Failure(failure) => handleIoFailure(command, failure)
    }
  }

  private def handleIoFailure(command: IoActorCommand[_], failure: Throwable) = {
    command.timeBeforeNextTry match {
      case None => command.fail(failure)
      case Some(interval) => handleFailure(command, failure, interval)
    }
  }

  private def handleFailure(command: IoActorCommand[_], failure: Throwable, retryInterval: FiniteDuration) = {
    if (isRetryable(failure)) {
      context.system.scheduler.scheduleOnce(retryInterval, context.parent, command.withNextBackoff)
      command.retried(failure)
    } else command.fail(failure)
  }
  
  private def isRetryable(throwable: Throwable) = throwable match {
      // Retries on 408, 429, 500, 502, 503, 504 and "internalError"
    case e: StorageException => e.retryable()
    case _ => false
  }

  /* Command implementations */
  private def copy(source: Path, destination: Path, overwrite: Boolean): Try[Unit] = Try {
    val destinationFile = File(destination)
    // Create directories
    if (destinationFile.isDirectory) destinationFile.createDirectories()
    else destinationFile.parentOption match {
      case Some(parent) => parent.createDirectories()
      case _ =>
    }
    
    File(source).copyTo(destination, overwrite = overwrite)
  }.void
  private def readAsString(file: Path): Try[String] = Try(File(file).contentAsString)
  private def size(file: Path): Try[Long] = Try(File(file).size)
  private def write(file: Path, content: String)(implicit openOptions: OpenOptions): Try[Unit] = Try(File(file).write(content)).void
  private def delete(file: Path, swallowIOExceptions: Boolean): Try[Unit] = Try(File(file).delete(swallowIOExceptions = swallowIOExceptions)).void
}

object NioServiceWorker {
  def props() = Props(new NioServiceWorker())
} 
