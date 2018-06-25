package cromwell.core.io

import akka.actor.ActorRef
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.io.IoPromiseProxyActor.IoCommandWithPromise
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.Path

import scala.concurrent.Future
import scala.concurrent.duration._
import net.ceedubs.ficus.Ficus._

object AsyncIo {
  private val ioTimeouts = ConfigFactory.load().as[Config]("system.io.timeout")
  val defaultTimeout = ioTimeouts.as[FiniteDuration]("default")
  val copyTimeout = ioTimeouts.as[FiniteDuration]("copy")
}

/**
  * Provides Futurized methods for I/O actions processed through the IoActor
  */
class AsyncIo(ioEndpoint: ActorRef, ioCommandBuilder: IoCommandBuilder) {
  private def asyncCommand[A](command: IoCommand[A], timeout: FiniteDuration = AsyncIo.defaultTimeout) = {
    val commandWithPromise = IoCommandWithPromise(command, timeout)
    ioEndpoint ! commandWithPromise
    commandWithPromise.promise.future
  }
  
  /**
    * IMPORTANT: This loads the entire content of the file into memory !
    * Only use for small files !
    */
  def contentAsStringAsync(path: Path, maxBytes: Option[Int], failOnOverflow: Boolean): Future[String] = {
    asyncCommand(ioCommandBuilder.contentAsStringCommand(path, maxBytes, failOnOverflow))
  }

  def writeAsync(path: Path, content: String, options: OpenOptions): Future[Unit] = {
    asyncCommand(ioCommandBuilder.writeCommand(path, content, options))
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

  def copyAsync(src: Path, dest: Path, overwrite: Boolean = true): Future[Unit] = {
    // Allow for a much larger timeout for copies, as large files can take a while (even on gcs, if they are in different locations...)
    asyncCommand(ioCommandBuilder.copyCommand(src, dest, overwrite), AsyncIo.copyTimeout)
  }
}
