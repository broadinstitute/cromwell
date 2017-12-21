package cromwell.core.io

import cromwell.core.io.DefaultIoCommand._
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.Path

/**
  * Can be used to customize IoCommands for the desired I/O operations
  */
abstract class PartialIoCommandBuilder {
  def contentAsStringCommand: PartialFunction[Path, IoContentAsStringCommand] = PartialFunction.empty
  def writeCommand: PartialFunction[(Path, String, OpenOptions), IoWriteCommand] = PartialFunction.empty
  def sizeCommand: PartialFunction[Path, IoSizeCommand] = PartialFunction.empty
  def deleteCommand: PartialFunction[(Path, Boolean), IoDeleteCommand] = PartialFunction.empty
  def copyCommand: PartialFunction[(Path, Path, Boolean), IoCopyCommand] = PartialFunction.empty
  def hashCommand: PartialFunction[Path, IoHashCommand] = PartialFunction.empty
  def touchCommand: PartialFunction[Path, IoTouchCommand] = PartialFunction.empty
  def existsCommand: PartialFunction[Path, IoExistsCommand] = PartialFunction.empty
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
  def contentAsStringCommand(path: Path): IoContentAsStringCommand = {
    partialBuilders.toStream.map(_.contentAsStringCommand.lift(path)).collectFirst({
      case Some(command) => command
    }).getOrElse(DefaultIoContentAsStringCommand(path))
  }
  
  def writeCommand(path: Path, content: String, options: OpenOptions): IoWriteCommand = {
    partialBuilders.toStream.map(_.writeCommand.lift((path, content, options))).collectFirst({
      case Some(command) => command
    }).getOrElse(DefaultIoWriteCommand(path, content, options))
  }
  
  def sizeCommand(path: Path): IoSizeCommand = {
    partialBuilders.toStream.map(_.sizeCommand.lift(path)).collectFirst({
      case Some(command) => command
    }).getOrElse(DefaultIoSizeCommand(path))
  } 
  
  def deleteCommand(path: Path, swallowIoExceptions: Boolean = true): IoDeleteCommand = {
    partialBuilders.toStream.map(_.deleteCommand.lift((path, swallowIoExceptions))).collectFirst({
      case Some(command) => command
    }).getOrElse(DefaultIoDeleteCommand(path, swallowIoExceptions))
  }
  
  def copyCommand(src: Path, dest: Path, overwrite: Boolean): IoCopyCommand = {
    partialBuilders.toStream.map(_.copyCommand.lift((src, dest, overwrite))).collectFirst({
      case Some(command) => command
    }).getOrElse(DefaultIoCopyCommand(src, dest, overwrite))
  }
  
  def hashCommand(file: Path): IoHashCommand = {
    partialBuilders.toStream.map(_.hashCommand.lift(file)).collectFirst({
      case Some(command) => command
    }).getOrElse(DefaultIoHashCommand(file))
  }
  
  def touchCommand(file: Path): IoTouchCommand = {
    partialBuilders.toStream.map(_.touchCommand.lift(file)).collectFirst({
      case Some(command) => command
    }).getOrElse(DefaultIoTouchCommand(file))
  }

  def existsCommand(file: Path): IoExistsCommand = {
    partialBuilders.toStream.map(_.existsCommand.lift(file)).collectFirst({
      case Some(command) => command
    }).getOrElse(DefaultIoExistsCommand(file))
  }

  def readLines(file: Path): IoReadLinesCommand = {
    partialBuilders.toStream.map(_.readLinesCommand.lift(file)).collectFirst({
      case Some(command) => command
    }).getOrElse(DefaultIoReadLinesCommand(file))
  }
}

/**
  * Only builds DefaultIoCommands.
  */
case object DefaultIoCommandBuilder extends IoCommandBuilder(List.empty)
