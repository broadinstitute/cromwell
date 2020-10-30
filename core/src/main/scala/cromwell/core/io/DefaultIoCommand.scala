package cromwell.core.io

import better.files.File.OpenOptions
import cromwell.core.io.IoContentAsStringCommand.IoReadOptions
import cromwell.core.path.Path

object DefaultIoCommand {
  case class DefaultIoCopyCommand(override val source: Path,
                                  override val destination: Path,
                                  override val overwrite: Boolean) extends IoCopyCommand(source, destination, overwrite) {
    override def commandDescription: String = s"DefaultIoCopyCommand source '$source' destination '$destination' overwrite '$overwrite'"
  }

  case class DefaultIoContentAsStringCommand(override val file: Path, override val options: IoReadOptions) extends IoContentAsStringCommand(file, options) {
    override def commandDescription: String = s"DefaultIoContentAsStringCommand file '$file' options '$options'"
  }

  case class DefaultIoSizeCommand(override val file: Path) extends IoSizeCommand(file) {
    override def commandDescription: String = s"DefaultIoSizeCommand file '$file'"
  }

  case class DefaultIoWriteCommand(override val file: Path,
                                   override val content: String,
                                   override val openOptions: OpenOptions,
                                   override val compressPayload: Boolean) extends IoWriteCommand(
    file, content, openOptions, compressPayload
  ) {
    override def commandDescription: String = s"DefaultIoWriteCommand file '$file' content length " +
      s"'${content.length}' openOptions '$openOptions' compressPayload '$compressPayload'"
  }

  case class DefaultIoDeleteCommand(override val file: Path,
                                    override val swallowIOExceptions: Boolean) extends IoDeleteCommand(
    file, swallowIOExceptions
  ) {
    override def commandDescription: String = s"DefaultIoDeleteCommand file '$file' swallowIOExceptions '$swallowIOExceptions'"
  }

  case class DefaultIoHashCommand(override val file: Path) extends IoHashCommand(file) {
    override def commandDescription: String = s"DefaultIoHashCommand file '$file'"
  }

  case class DefaultIoTouchCommand(override val file: Path) extends IoTouchCommand(file) {
    override def commandDescription: String = s"DefaultIoTouchCommand file '$file'"
  }

  case class DefaultIoExistsCommand(override val file: Path) extends IoExistsCommand(file) {
    override def commandDescription: String = s"DefaultIoExistsCommand file '$file'"
  }

  case class DefaultIoReadLinesCommand(override val file: Path) extends IoReadLinesCommand(file) {
    override def commandDescription: String = s"DefaultIoReadLinesCommand file '$file'"
  }

  case class DefaultIoIsDirectoryCommand(override val file: Path) extends IoIsDirectoryCommand(file) {
    override def commandDescription: String = s"DefaultIoIsDirectoryCommand file '$file'"
  }
}
