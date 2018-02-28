package cromwell.engine.instrumentation

import akka.actor.Actor
import cats.data.NonEmptyList
import cromwell.core.instrumentation.InstrumentationKeys._
import cromwell.core.instrumentation.InstrumentationPrefixes._
import cromwell.core.io._
import cromwell.engine.io.IoActor.IoResult
import cromwell.filesystems.gcs.{GcsPath, GoogleUtil}
import cromwell.services.instrumentation.CromwellInstrumentation._
import cromwell.services.instrumentation.CromwellInstrumentationActor

/**
  * Implicit methods to help convert Io objects to metric values
  */
private object IoInstrumentationImplicits {
  val LocalPath = NonEmptyList.of("local")
  val GcsPath = NonEmptyList.of("gcs")
  val UnknownFileSystemPath = NonEmptyList.of("unknown")
  
  val backpressure = NonEmptyList.of("backpressure")

  /**
    * Augments IoResult to provide instrumentation conversion methods
    */
  implicit class InstrumentedIoResult(val ioResult: IoResult) extends AnyVal {
    /**
      * Returns the instrumentation path of this IoResult
      */
    def toPath: InstrumentationPath = ioResult match {
      case (_: IoSuccess[_], ioCommandContext) => ioCommandContext.request.successPath
      case (f: IoFailure[_], ioCommandContext) => ioCommandContext.request.failedPath(f.failure)
    }
  }

  /**
    * Augments IoCommand to provide instrumentation conversion methods
    */
  implicit class InstrumentedIoCommand(val ioCommand: IoCommand[_]) extends AnyVal {
    /**
      * Returns the instrumentation path of this IoCommand
      */
    def toPath: InstrumentationPath = {
      val path = ioCommand match {
        case copy: IoCopyCommand => (copy.source, copy.destination) match {
          case (_: GcsPath, _) | (_, _: GcsPath) => GcsPath
          case _ => LocalPath
        }
        case singleFileCommand: SingleFileIoCommand[_] => singleFileCommand.file match {
          case _: GcsPath => GcsPath
          case _ => LocalPath
        }
        case _ => UnknownFileSystemPath
      }

      path.concat(ioCommand.name)
    }

    /**
      * Returns a successful instrumentation path for this IoCommand
      */
    def successPath: InstrumentationPath = ioCommand.toPath.concat(SuccessKey)

    /**
      * Returns a failed instrumentation path for this IoCommand provided a throwable
      */
    def failedPath(failure: Throwable): InstrumentationPath = {
      ioCommand.toPath.concat(FailureKey).withStatusCodeFailure(GoogleUtil.extractStatusCode(failure))
    }

    /**
      * Returns a retried instrumentation path for this IoCommand provided a throwable
      */
    def retriedPath(failure: Throwable): InstrumentationPath = {
      ioCommand.toPath.concat(RetryKey).withStatusCodeFailure(GoogleUtil.extractStatusCode(failure))
    }
  }
}

/**
  * Helper methods for Io instrumentation
  */
trait IoInstrumentation extends CromwellInstrumentationActor { this: Actor =>
  import IoInstrumentationImplicits._

  /**
    * Base increment method for IoInstrumentation
    */
  final private def incrementIo(statsDPath: InstrumentationPath): Unit = increment(statsDPath, IoPrefix)

  /**
    * Increment an IoResult to the proper bucket depending on the request type and the result (success or failure).
    */
  final def incrementIoResult(ioResult: IoResult): Unit = incrementIo(ioResult.toPath)
  
  final def incrementBackpressure(): Unit = incrementIo(backpressure)

  /**
    * Increment an IoCommand to the proper bucket depending on the request type.
    */
  final def incrementIoRetry(ioCommand: IoCommand[_], failure: Throwable): Unit = incrementIo(ioCommand.retriedPath(failure))
}
