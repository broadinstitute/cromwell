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
  */
class IoCopyCommand(val source: Path, val destination: Path, val overwrite: Boolean) extends IoCommand[Unit]

/**
  * Read file as a string (load the entire content in memory)
  */
class IoContentAsStringCommand(val file: Path) extends IoCommand[String]

/**
  * Return the size of file
  */
class IoSizeCommand(val file: Path) extends IoCommand[Long]

/**
  * Write content in file
  */
class IoWriteCommand(val file: Path, val content: String, val openOptions: OpenOptions) extends IoCommand[Unit]

/**
  * Delete file
  */
class IoDeleteCommand(val file: Path, val swallowIOExceptions: Boolean) extends IoCommand[Unit]

/**
  * Get Hash value for file
  */
class IoHashCommand(val file: Path) extends IoCommand[String]
