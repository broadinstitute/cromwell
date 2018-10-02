package cloud.nio.impl.ftp

import java.io.IOException
import java.nio.file.attribute.FileAttribute
import java.nio.file.{FileAlreadyExistsException, NoSuchFileException, Path}

import cloud.nio.impl.ftp.FtpUtil.FtpIoException
import cloud.nio.spi.{CloudNioFileSystemProvider, CloudNioPath, CloudNioReadChannel, CloudNioRetry}
import com.typesafe.config.Config
import com.typesafe.scalalogging.StrictLogging

import scala.util.{Failure, Success, Try}

class FtpCloudNioFileSystemProvider(override val config: Config, val credentials: FtpCredentials, ftpFileSystems: FtpFileSystems) extends CloudNioFileSystemProvider with StrictLogging {
  val ftpConfig = ftpFileSystems.config

  override def fileProvider = new FtpCloudNioFileProvider(this)

  override def isFatal(exception: Exception) = exception match {
    case ftpException: FtpIoException => ftpException.isFatal
    case _ => true
  }
  override def isTransient(exception: Exception) = false
  override def getScheme = "ftp"

  override def cloudNioReadChannel(retry: CloudNioRetry, cloudNioPath: CloudNioPath) = {
    /*
     * This is important, we need to get the file size before and give it to the read channel. Otherwise the read channel
     * will try to get it using the fileProvider which will require a new client lease and can result in a  deadlock of the client pool, since
     * the read channel holds on to its lease until its closed.
     */
    val preComputedFileSize = retry.from(() => fileProvider.fileAttributes(cloudNioPath.cloudHost, cloudNioPath.cloudPath).map(_.size()))
    new CloudNioReadChannel(fileProvider, retry, cloudNioPath) {
      override def fileSize = preComputedFileSize
    }
  }

  override def createDirectory(dir: Path, attrs: FileAttribute[_]*): Unit = {
    Try {
      retry.from(() => {
        val cloudNioPath = CloudNioPath.checkPath(dir)
        fileProvider.createDirectory(cloudNioPath.cloudHost, cloudNioPath.cloudPath)
      })
    } match {
      case Success(_) =>
      case Failure(f: FileAlreadyExistsException) => throw f
      case Failure(f: NoSuchFileException) => throw f
      // java.nio.files.Files.createDirectories swallows IOException, that are not FileAlreadyExistsException, so log them here so we can now what happened
      case Failure(f: IOException) =>
        logger.debug(s"Failed to create directory at $dir", f)
        throw f
      case Failure(f) => throw f
    }
  }

  override def usePseudoDirectories = false

  override def newCloudNioFileSystem(uriAsString: String, config: Config) = {
    newCloudNioFileSystemFromHost(getHost(uriAsString))
  }
  
  def newCloudNioFileSystemFromHost(host: String) = {
    ftpFileSystems.getFileSystem(host, this)
  }
}
