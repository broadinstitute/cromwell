package cromwell.engine.io.nio

import akka.actor.ActorSystem
import akka.stream.scaladsl.Flow
import cats.effect._
import com.typesafe.config.Config
import common.util.IORetry
import cromwell.core.io._
import cromwell.core.path.Path
import cromwell.engine.io.IoActor._
import cromwell.engine.io.RetryableRequestSupport.{isInfinitelyRetryable, isRetryable}
import cromwell.engine.io.{IoAttempts, IoCommandContext, IoCommandStalenessBackpressuring}
import cromwell.filesystems.http.HttpPath
import net.ceedubs.ficus.Ficus._
import net.ceedubs.ficus.readers.ValueReader

import java.io._
import java.nio.charset.StandardCharsets
import scala.concurrent.ExecutionContext
import scala.concurrent.duration.FiniteDuration

/**
  * Flow that executes IO operations by calling java.nio.Path methods
  */
class NioFlow(parallelism: Int,
              onRetryCallback: IoCommandContext[_] => Throwable => Unit,
              onBackpressure: Option[Double] => Unit,
              numberOfAttempts: Int,
              commandBackpressureStaleness: FiniteDuration
)(implicit system: ActorSystem)
    extends IoCommandStalenessBackpressuring {

  implicit private val ec: ExecutionContext = system.dispatcher
  implicit private val timer: Timer[IO] = IO.timer(ec)
  implicit private val contextShift: ContextShift[IO] = IO.contextShift(ec)

  override def maxStaleness: FiniteDuration = commandBackpressureStaleness

  private val processCommand: DefaultCommandContext[_] => IO[IoResult] = commandContext => {

    val onRetry: (Throwable, IoAttempts) => IoAttempts = (t, s) => {
      onRetryCallback(commandContext)
      IoAttempts.updateState(t, s)
    }

    val operationResult = IORetry.withRetry(
      handleSingleCommand(commandContext.request),
      IoAttempts(1),
      maxRetries = Option(numberOfAttempts),
      backoff = IoCommand.defaultBackoff,
      isRetryable = isRetryable,
      isInfinitelyRetryable = isInfinitelyRetryable,
      onRetry = onRetry
    )

    val io = for {
      _ <- IO.shift(ec)
      result <- operationResult
    } yield (result, commandContext)

    io handleErrorWith { failure =>
      IO.pure(commandContext.fail(failure))
    }
  }

  private[nio] def handleSingleCommand(ioSingleCommand: IoCommand[_]): IO[IoSuccess[_]] = {
    val ret = ioSingleCommand match {
      case copyCommand: IoCopyCommand => copy(copyCommand) map copyCommand.success
      case writeCommand: IoWriteCommand => write(writeCommand) map writeCommand.success
      case deleteCommand: IoDeleteCommand => delete(deleteCommand) map deleteCommand.success
      case sizeCommand: IoSizeCommand => size(sizeCommand) map sizeCommand.success
      case readAsStringCommand: IoContentAsStringCommand =>
        readAsString(readAsStringCommand) map readAsStringCommand.success
      case hashCommand: IoHashCommand => hash(hashCommand) map hashCommand.success
      case touchCommand: IoTouchCommand => touch(touchCommand) map touchCommand.success
      case existsCommand: IoExistsCommand => exists(existsCommand) map existsCommand.success
      case readLinesCommand: IoReadLinesCommand => readLines(readLinesCommand) map readLinesCommand.success
      case isDirectoryCommand: IoIsDirectoryCommand => isDirectory(isDirectoryCommand) map isDirectoryCommand.success
      case _ => IO.raiseError(new UnsupportedOperationException("Method not implemented"))
    }

    // Apply backpressure if this command has been waiting too long to execute.
    backpressureIfStale(ioSingleCommand, onBackpressure)
    ret
  }

  private[io] val flow =
    Flow[DefaultCommandContext[_]].mapAsyncUnordered[IoResult](parallelism)(processCommand.andThen(_.unsafeToFuture()))

  private def copy(copy: IoCopyCommand) = IO {
    createDirectories(copy.destination)
    copy.source.copyTo(copy.destination, overwrite = true)
    ()
  }

  private def write(write: IoWriteCommand) = IO {
    createDirectories(write.file)
    write.file.writeContent(write.content)(write.openOptions, StandardCharsets.UTF_8, write.compressPayload)
    ()
  }

  private def delete(delete: IoDeleteCommand) = IO {
    delete.file.delete(delete.swallowIOExceptions)
    ()
  }

  private def readAsString(command: IoContentAsStringCommand): IO[String] = IO {
    new String(
      command.file.limitFileContent(command.options.maxBytes, command.options.failOnOverflow),
      StandardCharsets.UTF_8
    ).replaceAll("\\r\\n", "\\\n")
  }

  private def size(size: IoSizeCommand) =
    size.file match {
      case httpPath: HttpPath => IO.fromFuture(IO(httpPath.fetchSize))
      case nioPath => IO(nioPath.size)
    }

  private def hash(hash: IoHashCommand): IO[String] =
    NioHashing.hash(hash.file, hash.hashStrategy)

  private def touch(touch: IoTouchCommand) = IO {
    touch.file.touch()
  }

  private def exists(exists: IoExistsCommand) = IO {
    exists.file.exists
  }

  private def readLines(exists: IoReadLinesCommand) = IO {
    exists.file.withReader { reader =>
      LazyList.continually(reader.readLine()).takeWhile(_ != null).toList
    }
  }

  private def isDirectory(isDirectory: IoIsDirectoryCommand) = IO {
    isDirectory.file.isDirectory
  }

  private def createDirectories(path: Path) = path.parent.createDirectories()
}

object NioFlow {
  case class NioFlowConfig(parallelism: Int)

  implicit val nioFlowConfigReader: ValueReader[NioFlowConfig] = (config: Config, path: String) => {
    val base = config.as[Config](path)
    val parallelism = base.as[Int]("parallelism")
    NioFlowConfig(parallelism)
  }
}

case class ChecksumFailedException(message: String) extends IOException(message)
