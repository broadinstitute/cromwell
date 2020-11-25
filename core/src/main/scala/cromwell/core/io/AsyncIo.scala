package cromwell.core.io

import akka.actor.ActorRef
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.Path

import scala.concurrent.Future
import scala.concurrent.duration._
import net.ceedubs.ficus.Ficus._

import scala.util.{Failure, Success, Try}

object AsyncIo {
  private val ioTimeouts = ConfigFactory.load().as[Config]("system.io.timeout")
  val defaultTimeout: FiniteDuration = ioTimeouts.as[FiniteDuration]("default")
  val copyTimeout: FiniteDuration = ioTimeouts.as[FiniteDuration]("copy")
}

/**
  * Provides Futurized methods for I/O actions processed through the IoActor
  */
class AsyncIo(ioEndpoint: ActorRef, ioCommandBuilder: IoCommandBuilder) {
  private def asyncCommand[A](commandTry: Try[IoCommand[A]],
                              timeout: FiniteDuration = AsyncIo.defaultTimeout): Future[A] = {
    commandTry match {
      case Failure(throwable) =>
        Future.failed(throwable)
      case Success(command) =>
        val commandWithPromise = IoCommandWithPromise(command, timeout)
        ioEndpoint ! commandWithPromise
        commandWithPromise.promise.future
    }
  }
  
  /**
    * IMPORTANT: This loads the entire content of the file into memory !
    * Only use for small files !
    */
  def contentAsStringAsync(path: Path, maxBytes: Option[Int], failOnOverflow: Boolean): Future[String] = {
    asyncCommand(ioCommandBuilder.contentAsStringCommand(path, maxBytes, failOnOverflow))
  }

  def writeAsync(path: Path, content: String, options: OpenOptions, compressPayload: Boolean = false): Future[Unit] = {
    asyncCommand(ioCommandBuilder.writeCommand(path, content, options, compressPayload))
  }

  def sizeAsync(path: Path): Future[Long] = {
    asyncCommand(ioCommandBuilder.sizeCommand(path))
  }

  def hashAsync(path: Path): Future[String] = {
    asyncCommand(ioCommandBuilder.hashCommand(path))
  }

  def deleteAsync(path: Path, swallowIoExceptions: Boolean = false): Future[Unit] = {
    asyncCommand(ioCommandBuilder.deleteCommand(path, swallowIoExceptions))
  }

  def existsAsync(path: Path): Future[Boolean] = {
    asyncCommand(ioCommandBuilder.existsCommand(path))
  }

  def readLinesAsync(path: Path): Future[Traversable[String]] = {
    asyncCommand(ioCommandBuilder.readLines(path))
  }
  
  def isDirectory(path: Path): Future[Boolean] = {
    asyncCommand(ioCommandBuilder.isDirectoryCommand(path))
  }

  def copyAsync(src: Path, dest: Path): Future[Unit] = {
    // Allow for a much larger timeout for copies, as large files can take a while (even on gcs, if they are in different locations...)
    asyncCommand(ioCommandBuilder.copyCommand(src, dest), AsyncIo.copyTimeout)
  }
}
