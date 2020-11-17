package cromwell.core.io

import cromwell.core.io.DefaultIoCommand._
import cromwell.core.io.IoContentAsStringCommand.IoReadOptions
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.Path

import scala.util.Try

/**
  * Can be used to customize IoCommands for the desired I/O operations
  */
//noinspection MutatorLikeMethodIsParameterless
abstract class PartialIoCommandBuilder {
  def contentAsStringCommand: PartialFunction[(Path, Option[Int], Boolean), Try[IoContentAsStringCommand]] =
    PartialFunction.empty
  def writeCommand: PartialFunction[(Path, String, OpenOptions, Boolean), Try[IoWriteCommand]] = PartialFunction.empty
  def sizeCommand: PartialFunction[Path, Try[IoSizeCommand]] = PartialFunction.empty
  def deleteCommand: PartialFunction[(Path, Boolean), Try[IoDeleteCommand]] = PartialFunction.empty
  def copyCommand: PartialFunction[(Path, Path), Try[IoCopyCommand]] = PartialFunction.empty
  def hashCommand: PartialFunction[Path, Try[IoHashCommand]] = PartialFunction.empty
  def touchCommand: PartialFunction[Path, Try[IoTouchCommand]] = PartialFunction.empty
  def existsCommand: PartialFunction[Path, Try[IoExistsCommand]] = PartialFunction.empty
  def isDirectoryCommand: PartialFunction[Path, Try[IoIsDirectoryCommand]] = PartialFunction.empty
  def readLinesCommand: PartialFunction[Path, Try[IoReadLinesCommand]] = PartialFunction.empty
}

object IoCommandBuilder {
  def apply(partialBuilders: PartialIoCommandBuilder*): IoCommandBuilder = {
    new IoCommandBuilder(partialBuilders.toList)
  }

  def apply: IoCommandBuilder = {
    new IoCommandBuilder(List.empty)
  }
}

/**
  * A command builder is a way to build an IoCommand for a specific I/O action.
  * The default IO command builder can execute all IoCommand and therefore can always be used.
  * One might want to create different I/O commands to allow for optimizations when the commands are processed by the I/O actor.
  * Currently the only other command builder is the GcsBatchCommandBuilder that overrides some of the operations
  * to return GcsBatchCommands instead that will be optimized by the IoActor.
  * 
  * This always defaults to building a DefaultIoCommand.
  * @param partialBuilders list of PartialIoCommandBuilder to try
  */
class IoCommandBuilder(partialBuilders: List[PartialIoCommandBuilder] = List.empty) {
  // Find the first partialBuilder for which the partial function is defined, or use the default
  private def buildOrDefault[A, B](builder: PartialIoCommandBuilder => PartialFunction[A, Try[B]],
                                   params: A,
                                   default: => B): Try[B] = {
    partialBuilders.toStream.map(builder(_).lift(params)).collectFirst({
      case Some(command) => command
    }).getOrElse(Try(default))
  }
  
  def contentAsStringCommand(path: Path,
                             maxBytes: Option[Int],
                             failOnOverflow: Boolean): Try[IoContentAsStringCommand] = {
    buildOrDefault(_.contentAsStringCommand, (path, maxBytes, failOnOverflow), DefaultIoContentAsStringCommand(path, IoReadOptions(maxBytes, failOnOverflow)))
  }
  
  def writeCommand(path: Path,
                   content: String,
                   options: OpenOptions,
                   compressPayload: Boolean = false): Try[IoWriteCommand] = {
    buildOrDefault(_.writeCommand, (path, content, options, compressPayload), DefaultIoWriteCommand(path, content, options, compressPayload))
  }
  
  def sizeCommand(path: Path): Try[IoSizeCommand] = {
    buildOrDefault(_.sizeCommand, path, DefaultIoSizeCommand(path))
  } 
  
  def deleteCommand(path: Path, swallowIoExceptions: Boolean = true): Try[IoDeleteCommand] = {
    buildOrDefault(_.deleteCommand, (path, swallowIoExceptions), DefaultIoDeleteCommand(path, swallowIoExceptions))
  }
  
  def copyCommand(src: Path, dest: Path): Try[IoCopyCommand] = {
    buildOrDefault(_.copyCommand, (src, dest), DefaultIoCopyCommand(src, dest))
  }
  
  def hashCommand(file: Path): Try[IoHashCommand] = {
    buildOrDefault(_.hashCommand, file, DefaultIoHashCommand(file))
  }
  
  def touchCommand(file: Path): Try[IoTouchCommand] = {
    buildOrDefault(_.touchCommand, file, DefaultIoTouchCommand(file))
  }

  def existsCommand(file: Path): Try[IoExistsCommand] = {
    buildOrDefault(_.existsCommand, file, DefaultIoExistsCommand(file))
  }

  def isDirectoryCommand(file: Path): Try[IoIsDirectoryCommand] = {
    buildOrDefault(_.isDirectoryCommand, file, DefaultIoIsDirectoryCommand(file))
  }

  def readLines(file: Path): Try[IoReadLinesCommand] = {
    buildOrDefault(_.readLinesCommand, file, DefaultIoReadLinesCommand(file))
  }
}

/**
  * Only builds DefaultIoCommands.
  */
case object DefaultIoCommandBuilder extends IoCommandBuilder(List.empty)
