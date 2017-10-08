package cromwell.core.io

import cromwell.core.io.DefaultIoCommand._
import cromwell.core.path.BetterFileMethods.OpenOptions
import cromwell.core.path.Path

trait IoCommandBuilder {
  def contentAsStringCommand(path: Path): IoContentAsStringCommand
  def writeCommand(path: Path, content: String, options: OpenOptions): IoWriteCommand
  def sizeCommand(path: Path): IoSizeCommand
  def deleteCommand(path: Path, swallowIoExceptions: Boolean): IoDeleteCommand
  def copyCommand(src: Path, dest: Path, overwrite: Boolean): IoCopyCommand
  def hashCommand(file: Path): IoHashCommand
  def touchCommand(file: Path): IoTouchCommand
}

trait DefaultIoCommandBuilder extends IoCommandBuilder {
  def contentAsStringCommand(path: Path): IoContentAsStringCommand = DefaultIoContentAsStringCommand(path)
  def writeCommand(path: Path, content: String, options: OpenOptions): IoWriteCommand = DefaultIoWriteCommand(path, content, options)
  def sizeCommand(path: Path): IoSizeCommand = DefaultIoSizeCommand(path)
  def deleteCommand(path: Path, swallowIoExceptions: Boolean): IoDeleteCommand = DefaultIoDeleteCommand(path, swallowIoExceptions)
  def copyCommand(src: Path, dest: Path, overwrite: Boolean): IoCopyCommand = DefaultIoCopyCommand(src, dest, overwrite)
  def hashCommand(file: Path): IoHashCommand = DefaultIoHashCommand(file)
  def touchCommand(file: Path): IoTouchCommand = DefaultIoTouchCommand(file)
}
