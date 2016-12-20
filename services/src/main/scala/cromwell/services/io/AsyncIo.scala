package cromwell.services.io

import java.nio.file.Path

import akka.actor.ActorRef
import better.files.File.OpenOptions
import cromwell.services.io.promise._

import scala.concurrent.{Future, Promise}

trait AsyncIo {
  def serviceRegistryActor: ActorRef

  def size(file: Path): Future[Long] = {
    val promise = Promise[Long]
    serviceRegistryActor ! SizeCommandPromise(file)(promise)
    promise.future
  }
  
  def read(file: Path): Future[String] = {
    val promise = Promise[String]
    serviceRegistryActor ! ReadAsStringCommandPromise(file)(promise)
    promise.future
  }

  def write(file: Path, content: String, openOptions: OpenOptions = OpenOptions.default): Future[Unit] = {
    val promise = Promise[Unit]
    serviceRegistryActor ! WriteCommandPromise(file, content, openOptions)(promise)
    promise.future
  }
  
  def copy(source: Path, destination: Path, overwrite: Boolean = true): Future[Unit] = {
    val promise = Promise[Unit]
    serviceRegistryActor ! CopyCommandPromise(source, destination, overwrite)(promise)
    promise.future
  }

  def delete(source: Path, swallowIOExceptions: Boolean = true): Future[Unit] = {
    val promise = Promise[Unit]
    serviceRegistryActor ! DeleteCommandPromise(source, swallowIOExceptions)(promise)
    promise.future
  }
}
