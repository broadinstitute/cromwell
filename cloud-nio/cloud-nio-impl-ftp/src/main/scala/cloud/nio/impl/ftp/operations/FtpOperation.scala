package cloud.nio.impl.ftp.operations

import java.io.{InputStream, OutputStream}
import java.nio.file.{FileAlreadyExistsException, NoSuchFileException}

import cats.effect.IO
import cloud.nio.impl.ftp.FtpUtil._
import cloud.nio.spi.CloudNioFileSystem
import org.apache.commons.net.ftp.{FTPClient, FTPFile, FTPReply}

sealed trait FtpOperation[A] {
  def cloudHost: String
  def cloudPath: String
  def description: String
  def action: FTPClient => A
  def run(client: FTPClient): IO[A]

  def fullPath = s"ftp://$cloudHost/$cloudPath"
  protected def errorMessage = s"Failed to $description at $fullPath"

  protected def fail(client: FTPClient, cause: Option[Throwable] = None) =
    IO.raiseError[A](generateException(client, cause))

  private[operations] def generateException(client: FTPClient, cause: Option[Throwable]) = cause match {
    case None
        if client.getReplyCode == FTPReply.FILE_UNAVAILABLE && client.getReplyString.toLowerCase.contains("exists") =>
      new FileAlreadyExistsException(fullPath)
    case None
        if client.getReplyCode == FTPReply.FILE_UNAVAILABLE && client.getReplyString.toLowerCase.contains(
          "no such file"
        ) =>
      new NoSuchFileException(fullPath)
    case None =>
      FtpIoException(errorMessage, client.getReplyCode, Option(client.getReplyString).getOrElse("N/A"))
    case Some(c) =>
      FtpIoException(errorMessage, client.getReplyCode, Option(client.getReplyString).getOrElse("N/A"), Option(c))
  }

  protected def handleError(client: FTPClient)(failure: Throwable) = fail(client, Option(failure))

  protected def commonRun(client: FTPClient, bind: A => IO[A]): IO[A] =
    IO(action(client)) redeemWith (handleError(client), bind)

  override def toString = s"$description at $fullPath"
}

sealed trait FtpBooleanOperation extends FtpOperation[Boolean] {
  protected def failOnFalse: Boolean = true

  def run(client: FTPClient): IO[Boolean] =
    commonRun(client,
              {
                // Operation didn't throw but the result is false which means it failed
                case false if failOnFalse => fail(client)
                case result => IO.pure(result)
              }
    )
}

sealed trait FtpValueOperation[A <: AnyRef] extends FtpOperation[A] {
  def run(client: FTPClient): IO[A] =
    commonRun(client,
              {
                // Operation didn't throw but the result is null which means it failed
                case null => fail(client)
                case result => IO.pure(result)
              }
    )
}

case class FtpListFiles(cloudHost: String, cloudPath: String, description: String = "List files")
    extends FtpValueOperation[Array[FTPFile]] {
  override val action = _.listFiles(cloudPath.ensureSlashedPrefix)
}

case class FtpListDirectories(cloudHost: String, cloudPath: String, description: String = "List files")
    extends FtpValueOperation[Array[FTPFile]] {
  // We need to list the directories in the parent and see if any matches the name, hence the string manipulations
  lazy val parts =
    cloudPath.ensureSlashedPrefix.stripSuffix(CloudNioFileSystem.Separator).split(CloudNioFileSystem.Separator)
  lazy val parent = parts.init.mkString(CloudNioFileSystem.Separator)
  lazy val directoryName = parts.last

  override val action = _.listDirectories(parent)
}

case class FtpDeleteFile(cloudHost: String, cloudPath: String, description: String = "delete file")
    extends FtpBooleanOperation {
  override val action = _.deleteFile(cloudPath.ensureSlashedPrefix)
  override val failOnFalse = false
}

case class FtpInputStream(cloudHost: String, cloudPath: String, offset: Long, description: String = "read")
    extends FtpValueOperation[InputStream] {
  override val action = { client =>
    client.setRestartOffset(offset)
    client.retrieveFileStream(cloudPath.ensureSlashedPrefix)
  }
}

case class FtpOutputStream(cloudHost: String, cloudPath: String, description: String = "write")
    extends FtpValueOperation[OutputStream] {
  override val action = _.storeFileStream(cloudPath.ensureSlashedPrefix)
}

case class FtpCreateDirectory(cloudHost: String, cloudPath: String, description: String = "create directory")
    extends FtpBooleanOperation {
  override val action = _.makeDirectory(cloudPath.ensureSlashedPrefix)
}

case class FtpCompletePendingCommand(cloudHost: String, cloudPath: String, description: String = "close stream")
    extends FtpBooleanOperation {
  override val action = _.completePendingCommand()
}
