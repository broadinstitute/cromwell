package cromwell.core.io

import better.files.File.OpenOptions
import com.google.api.client.util.ExponentialBackOff
import cromwell.core.path.Path
import cromwell.core.retry.SimpleExponentialBackoff

import scala.concurrent.duration.{FiniteDuration, _}
import scala.language.postfixOps

object IoCommand {
  def defaultGoogleBackoff = new ExponentialBackOff.Builder()
    .setInitialIntervalMillis((1 second).toMillis.toInt)
    .setMaxIntervalMillis((5 minutes).toMillis.toInt)
    .setMultiplier(3L)
    .setRandomizationFactor(0.2D)
    .setMaxElapsedTimeMillis((10 minutes).toMillis.toInt)
    .build()
  
  def defaultBackoff = SimpleExponentialBackoff(defaultGoogleBackoff)
  
  type RetryCommand[T] = (FiniteDuration, IoCommand[T])
}

trait IoCommand[+T] {
  /**
    * Completes the command successfully
    * @return a message to be sent back to the sender, if needed
    */
  def success[S >: T](value: S): IoSuccess[S] = IoSuccess(this, value)
  
  /**
    * Fail the command with an exception
    */
  def fail[S >: T](failure: Throwable): IoFailure[S] = IoFailure(this, failure)
}

/**
  * Copy source -> destination
  * Will create the destination directory if it doesn't exist.
  */
class IoCopyCommand(val source: Path, val destination: Path, val overwrite: Boolean) extends IoCommand[Unit] {
  override def toString = s"copy ${source.pathAsString} to ${destination.pathAsString} with overwrite = $overwrite"
}

/**
  * Read file as a string (load the entire content in memory)
  */
class IoContentAsStringCommand(val file: Path) extends IoCommand[String] {
  override def toString = s"read content of ${file.pathAsString}"
}

/**
  * Return the size of file
  */
class IoSizeCommand(val file: Path) extends IoCommand[Long] {
  override def toString = s"get size of ${file.pathAsString}"
}

/**
  * Write content in file
  * Will create the destination directory if it doesn't exist.
  */
class IoWriteCommand(val file: Path, val content: String, val openOptions: OpenOptions) extends IoCommand[Unit] {
  override def toString = s"write to ${file.pathAsString}"
}

/**
  * Delete file
  */
class IoDeleteCommand(val file: Path, val swallowIOExceptions: Boolean) extends IoCommand[Unit] {
  override def toString = s"delete ${file.pathAsString}"
}

/**
  * Get Hash value for file
  */
class IoHashCommand(val file: Path) extends IoCommand[String] {
  override def toString = s"get hash of ${file.pathAsString}"
}
