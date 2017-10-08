package cromwell.core.io

import better.files.File.OpenOptions
import cromwell.core.path.Path

object DefaultIoCommand {
  case class DefaultIoCopyCommand(override val source: Path,
                                  override val destination: Path,
                                  override val overwrite: Boolean) extends IoCopyCommand(
    source, destination, overwrite
  )
  case class DefaultIoContentAsStringCommand(override val file: Path) extends IoContentAsStringCommand(file)
  case class DefaultIoSizeCommand(override val file: Path) extends IoSizeCommand(file)
  case class DefaultIoWriteCommand(override val file: Path,
                                   override val content: String,
                                   override val openOptions: OpenOptions) extends IoWriteCommand(
    file, content, openOptions
  )
  case class DefaultIoDeleteCommand(override val file: Path,
                                    override val swallowIOExceptions: Boolean) extends IoDeleteCommand(
    file, swallowIOExceptions
  )
  case class DefaultIoHashCommand(override val file: Path) extends IoHashCommand(file)
  case class DefaultIoTouchCommand(override val file: Path) extends IoTouchCommand(file)
}
