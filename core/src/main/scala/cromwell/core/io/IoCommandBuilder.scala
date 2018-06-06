package cromwell.core.io

import cromwell.core.io.DefaultIoCommand._
import cromwell.core.io.IoContentAsStringCommand.IoReadOptions
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.Path

/**
  * Can be used to customize IoCommands for the desired I/O operations
  */
abstract class PartialIoCommandBuilder {
  def contentAsStringCommand: PartialFunction[(Path, Option[Int], Boolean), IoContentAsStringCommand] = PartialFunction.empty
  def writeCommand: PartialFunction[(Path, String, OpenOptions), IoWriteCommand] = PartialFunction.empty
  def sizeCommand: PartialFunction[Path, IoSizeCommand] = PartialFunction.empty
  def deleteCommand: PartialFunction[(Path, Boolean), IoDeleteCommand] = PartialFunction.empty
  def copyCommand: PartialFunction[(Path, Path, Boolean), IoCopyCommand] = PartialFunction.empty
  def hashCommand: PartialFunction[Path, IoHashCommand] = PartialFunction.empty
  def touchCommand: PartialFunction[Path, IoTouchCommand] = PartialFunction.empty
  def existsCommand: PartialFunction[Path, IoExistsCommand] = PartialFunction.empty
  def isDirectoryCommand: PartialFunction[Path, IoIsDirectoryCommand] = PartialFunction.empty
  def readLinesCommand: PartialFunction[Path, IoReadLinesCommand] = PartialFunction.empty
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
  private def buildOrDefault[A, B](builder: PartialIoCommandBuilder => PartialFunction[A, B],
                                   params: A,
                                   default: B) = {
    partialBuilders.toStream.map(builder(_).lift(params)).collectFirst({
      case Some(command) => command
    }).getOrElse(default)
  }
  
  def contentAsStringCommand(path: Path, maxBytes: Option[Int], failOnOverflow: Boolean): IoContentAsStringCommand = {
    buildOrDefault(_.contentAsStringCommand, (path, maxBytes, failOnOverflow), DefaultIoContentAsStringCommand(path, IoReadOptions(maxBytes, failOnOverflow)))
  }
  
  def writeCommand(path: Path, content: String, options: OpenOptions): IoWriteCommand = {
    buildOrDefault(_.writeCommand, (path, content, options), DefaultIoWriteCommand(path, content, options))
  }
  
  def sizeCommand(path: Path): IoSizeCommand = {
    buildOrDefault(_.sizeCommand, path, DefaultIoSizeCommand(path))
  } 
  
  def deleteCommand(path: Path, swallowIoExceptions: Boolean = true): IoDeleteCommand = {
    buildOrDefault(_.deleteCommand, (path, swallowIoExceptions), DefaultIoDeleteCommand(path, swallowIoExceptions))
  }
  
  def copyCommand(src: Path, dest: Path, overwrite: Boolean): IoCopyCommand = {
    buildOrDefault(_.copyCommand, (src, dest, overwrite), DefaultIoCopyCommand(src, dest, overwrite))
  }
  
  def hashCommand(file: Path): IoHashCommand = {
    buildOrDefault(_.hashCommand, file, DefaultIoHashCommand(file))
  }
  
  def touchCommand(file: Path): IoTouchCommand = {
    buildOrDefault(_.touchCommand, file, DefaultIoTouchCommand(file))
  }

  def existsCommand(file: Path): IoExistsCommand = {
    buildOrDefault(_.existsCommand, file, DefaultIoExistsCommand(file))
  }

  def isDirectoryCommand(file: Path): IoIsDirectoryCommand = {
    buildOrDefault(_.isDirectoryCommand, file, DefaultIoIsDirectoryCommand(file))
  }

  def readLines(file: Path): IoReadLinesCommand = {
    buildOrDefault(_.readLinesCommand, file, DefaultIoReadLinesCommand(file))
  }
}

/**
  * Only builds DefaultIoCommands.
  */
case object DefaultIoCommandBuilder extends IoCommandBuilder(List.empty)
