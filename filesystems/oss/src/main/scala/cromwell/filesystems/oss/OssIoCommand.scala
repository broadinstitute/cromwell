package cromwell.filesystems.oss

import java.nio.file.OpenOption

import better.files.File.OpenOptions
import cromwell.core.io._

case object UploadFileOption extends OpenOption
case object UploadStringOption extends OpenOption

trait OssIoCommand[T] extends IoCommand[T]

case class OssIoCopyCommand(override val source: OssPath,
                            override val destination: OssPath,
                            override val overwrite: Boolean) extends IoCopyCommand(source, destination, overwrite) with OssIoCommand[Unit]

case class OssIoContentAsStringCommand(override val file: OssPath) extends IoContentAsStringCommand(file) with OssIoCommand[String]
case class OssIoSizeCommand(override val file: OssPath) extends IoSizeCommand(file) with OssIoCommand[Long]
case class OssIoWriteCommand(override val file: OssPath,
                             override val content: String,
                             override val openOptions: OpenOptions) extends IoWriteCommand(file, content, openOptions) with OssIoCommand[Unit]

case class OssIoDeleteCommand(override val file: OssPath,
                              override val swallowIOExceptions: Boolean) extends IoDeleteCommand(file, swallowIOExceptions) with OssIoCommand[Unit]

case class OssIoHashCommand(override val file: OssPath) extends IoHashCommand(file) with OssIoCommand[String]
case class OssIoTouchCommand(override val file: OssPath) extends IoTouchCommand(file) with OssIoCommand[Unit]
