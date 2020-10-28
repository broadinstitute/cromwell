package cromwell.core.io

import better.files.File.OpenOptions
import cromwell.core.io.IoContentAsStringCommand.IoReadOptions
import cromwell.core.path.Path

object DefaultIoCommand {
  case class DefaultIoCopyCommand(override val source: Path,
                                  override val destination: Path,
                                  override val overwrite: Boolean) extends IoCopyCommand(
    source, destination, overwrite
  ) {
    customDebug(s"DefaultIoCopyCommand.init source '$source' destination '$destination' overwrite '$overwrite'")
  }
  case class DefaultIoContentAsStringCommand(override val file: Path, override val options: IoReadOptions) extends IoContentAsStringCommand(file, options) {
    customDebug(s"DefaultIoCopyCommand.init file '$file' options '$options'")
  }
  case class DefaultIoSizeCommand(override val file: Path) extends IoSizeCommand(file) {
    customDebug(s"DefaultIoSizeCommand.init file '$file'")
  }
  case class DefaultIoWriteCommand(override val file: Path,
                                   override val content: String,
                                   override val openOptions: OpenOptions,
                                   override val compressPayload: Boolean) extends IoWriteCommand(
    file, content, openOptions, compressPayload
  ) {
    customDebug(s"DefaultIoWriteCommand.init file '$file' content length '${content.length}' openOptions '$openOptions' compressPayload '$compressPayload'")
  }

  case class DefaultIoDeleteCommand(override val file: Path,
                                    override val swallowIOExceptions: Boolean) extends IoDeleteCommand(
    file, swallowIOExceptions
  ) {
    customDebug(s"DefaultIoDeleteCommand.init file '$file' swallowIOExceptions '$swallowIOExceptions'")
  }
  case class DefaultIoHashCommand(override val file: Path) extends IoHashCommand(file) {
    customDebug(s"DefaultIoHashCommand.init file '$file'")
  }
  case class DefaultIoTouchCommand(override val file: Path) extends IoTouchCommand(file) {
    customDebug(s"DefaultIoTouchCommand.init file '$file'")
  }
  case class DefaultIoExistsCommand(override val file: Path) extends IoExistsCommand(file) {
    customDebug(s"DefaultIoExistsCommand.init file '$file'")
  }
  case class DefaultIoReadLinesCommand(override val file: Path) extends IoReadLinesCommand(file) {
    customDebug(s"DefaultIoReadLinesCommand.init file '$file'")
  }
  case class DefaultIoIsDirectoryCommand(override val file: Path) extends IoIsDirectoryCommand(file) {
    customDebug(s"DefaultIoIsDirectoryCommand.init file '$file'")
  }
}
