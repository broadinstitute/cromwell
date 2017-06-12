package cromwell.core.io

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
  def contentAsStringCommand(path: Path) = new IoContentAsStringCommand(path)
  def writeCommand(path: Path, content: String, options: OpenOptions) = new IoWriteCommand(path, content, options)
  def sizeCommand(path: Path) = new IoSizeCommand(path)
  def deleteCommand(path: Path, swallowIoExceptions: Boolean) = new IoDeleteCommand(path, swallowIoExceptions)
  def copyCommand(src: Path, dest: Path, overwrite: Boolean) = new IoCopyCommand(src, dest, overwrite)
  def hashCommand(file: Path) = new IoHashCommand(file)
  def touchCommand(file: Path) = new IoTouchCommand(file)
}
