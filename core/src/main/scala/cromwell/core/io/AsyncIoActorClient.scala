package cromwell.core.io

import akka.actor.ActorRef
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.Path

import scala.concurrent.Future

/**
  * Helper trait to used async Io method easily
  */
trait AsyncIoActorClient {
  
  def ioActor: ActorRef
  def ioCommandBuilder: IoCommandBuilder
  private val asyncIo = new AsyncIo(ioActor, ioCommandBuilder)
  
  /**
    * IMPORTANT: This loads the entire content of the file into memory !
    * Only use for small files !
    */
  def contentAsStringAsync(path: Path): Future[String] = {
    asyncIo.contentAsStringAsync(path)
  }
  def writeAsync(path: Path, content: String, options: OpenOptions): Future[Unit] = {
    asyncIo.writeAsync(path, content, options)
  }
  def sizeAsync(path: Path): Future[Long] = {
    asyncIo.sizeAsync(path)
  }
  def hashAsync(path: Path): Future[String] = {
    asyncIo.hashAsync(path)
  }
  def deleteAsync(path: Path, swallowIoExceptions: Boolean = false): Future[Unit] = {
    asyncIo.deleteAsync(path, swallowIoExceptions)
  }
  def existsAsync(path: Path): Future[Boolean] = {
    asyncIo.existsAsync(path)
  }
  def copyAsync(src: Path, dest: Path, overwrite: Boolean = true): Future[Unit] = {
    asyncIo.copyAsync(src, dest, overwrite)
  }
}
