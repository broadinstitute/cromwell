package cloud.nio.impl.ftp

import java.nio.file.Path
import java.nio.file.attribute.FileAttribute

import cloud.nio.impl.ftp.FtpUtil.FtpIoException
import cloud.nio.spi.{CloudNioFileSystemProvider, CloudNioPath, CloudNioReadChannel, CloudNioRetry}
import com.typesafe.config.Config

class FtpCloudNioFileSystemProvider(override val config: Config, val credentials: FtpCredentials) extends CloudNioFileSystemProvider {
  override def fileProvider = new FtpCloudNioFileProvider(this)

  override def isFatal(exception: Exception) = exception match {
    case ftpException: FtpIoException => ftpException.isFatal
    case _ => false
  }
  override def isTransient(exception: Exception) = false
  override def getScheme = "ftp"

  override def cloudNioReadChannel(retry: CloudNioRetry, cloudNioPath: CloudNioPath) = {
    /*
     * This is important, we need to get the file size before and give it to the read channel. Otherwise the read channel
     * will try to get it using the fileProvider which will require a new client lease and can result in a  deadlock of the client pool, since
     * the read channel holds on to its lease until its closed.
     */
    val preComputedFileSize = fileProvider.fileAttributes(cloudNioPath.cloudHost, cloudNioPath.cloudPath).map(_.size())
    new CloudNioReadChannel(fileProvider, retry, cloudNioPath) {
      override def fileSize = preComputedFileSize
    }
  }

  override def createDirectory(dir: Path, attrs: FileAttribute[_]*): Unit = {
    val cloudNioPath = CloudNioPath.checkPath(dir)
    fileProvider.createDirectory(cloudNioPath.cloudHost, cloudNioPath.cloudPath)
  }

  override def usePseudoDirectories = false
}
