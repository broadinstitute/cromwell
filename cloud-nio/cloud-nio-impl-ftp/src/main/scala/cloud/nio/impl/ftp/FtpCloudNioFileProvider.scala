package cloud.nio.impl.ftp

import java.io.{FileNotFoundException, IOException, InputStream, OutputStream}
import java.nio.channels.Channels
import java.nio.file.attribute.FileTime

import cloud.nio.impl.ftp.FtpCloudNioFileProvider._
import cloud.nio.impl.ftp.FtpFileSystems._
import cloud.nio.impl.ftp.FtpUtil._
import cloud.nio.spi.{CloudNioFileList, CloudNioFileProvider, CloudNioRegularFileAttributes}
import cloud.nio.util.TryWithResource._
import io.github.andrebeat.pool.Lease
import org.apache.commons.net.ftp.FTPClient
import org.apache.commons.net.io.Util

import scala.util.{Failure, Success}

object FtpCloudNioFileProvider {
  implicit class EnhancedCloudPath(val cloudPath: String) extends AnyVal {
    def ensureSlashPrefix = if (cloudPath.startsWith("/")) cloudPath else s"/$cloudPath"
  }
}

class FtpCloudNioFileProvider(fsProvider: FtpCloudNioFileSystemProvider) extends CloudNioFileProvider {
  override def existsPath(cloudHost: String, cloudPath: String) = {
    run(FtpOperation(cloudHost, cloudPath, "determine file existence"), acquireLease(cloudHost))(_.listFiles(cloudPath.ensureSlashPrefix)).nonEmpty
  }

  override def existsPaths(cloudHost: String, cloudPathPrefix: String) = {
    val maybeDirectory = run(FtpOperation(cloudHost, cloudPathPrefix, "determine files existence from prefix"), acquireLease(cloudHost))(_.listFiles(cloudPathPrefix.ensureSlashPrefix)).headOption
    maybeDirectory.exists(_.isDirectory)
  }

  /**
    * Returns a listing of keys within a bucket starting with prefix. The returned keys should include the prefix. The
    * paths must be absolute, but the key should not begin with a slash.
    */
  override def listObjects(cloudHost: String, cloudPathPrefix: String, markerOption: Option[String]) = {
    val files = run(FtpOperation(cloudHost, cloudPathPrefix, "list objects"), acquireLease(cloudHost)) { client =>
      client.listFiles(cloudPathPrefix.ensureSlashPrefix)
    }

    CloudNioFileList(files.map(_.getName), markerOption)
  }

  override def copy(sourceCloudHost: String, sourceCloudPath: String, targetCloudHost: String, targetCloudPath: String) = {
    val fileSystem = findFileSystem(sourceCloudHost)

    val (inputLease, outputLease) = fileSystem.leaseClientsPair

    val streams = () => {
      val is = inputStream(sourceCloudHost, sourceCloudPath, 0, inputLease)
      val os = outputStream(targetCloudHost, targetCloudPath, outputLease)
      new InputOutputStreams(is, os)
    }

    tryWithResource(streams) { ios =>
      Util.copyStream(ios.inputStream, ios.outputStream)
    } match {
      case Success(_) =>
      case Failure(failure) => throw new IOException(s"Failed to copy ftp://$sourceCloudHost/$sourceCloudPath to ftp://$targetCloudHost/$targetCloudPath", failure)
    }
  }

  override def deleteIfExists(cloudHost: String, cloudPath: String) = {
    runBoolean(FtpOperation(cloudHost, cloudPath, "delete"), acquireLease(cloudHost))(_.deleteFile(cloudPath.ensureSlashPrefix))
    true
  }

  private def inputStream(cloudHost: String, cloudPath: String, offset: Long, lease: Lease[FTPClient]) = {
    val operation = FtpOperation(cloudHost, cloudPath, "read")

    val stream: InputStream = run(operation, lease, autoRelease = false) { client =>
      client.setRestartOffset(offset)
      client.retrieveFileStream(cloudPath.ensureSlashPrefix)
    }

    // Wrap the input stream in a LeasedInputStream so that the lease can be release when the stream is closed
    new LeasedInputStream(cloudHost, cloudPath, stream, lease)
  }

  override def read(cloudHost: String, cloudPath: String, offset: Long) = {
    Channels.newChannel(inputStream(cloudHost, cloudPath, offset, acquireLease(cloudHost)))
  }

  private def outputStream(cloudHost: String, cloudPath: String, lease: Lease[FTPClient]) = {
    val operation = FtpOperation(cloudHost, cloudPath, "write")

    val stream: OutputStream = run(operation, lease, autoRelease = false)(_.storeFileStream(cloudPath.ensureSlashPrefix))

    // Wrap the input stream in a LeasedOutputStream so that the lease can be release when the stream is closed
    new LeasedOutputStream(cloudHost: String, cloudPath: String, stream, lease)
  }

  override def write(cloudHost: String, cloudPath: String) = {
    Channels.newChannel(outputStream(cloudHost: String, cloudPath: String, acquireLease(cloudHost)))
  }

  override def fileAttributes(cloudHost: String, cloudPath: String) = {
    val file = run(FtpOperation(cloudHost, cloudPath, "get file attributes"), acquireLease(cloudHost))(_.listFiles(cloudPath.ensureSlashPrefix)).headOption

    Option {
      new CloudNioRegularFileAttributes {
        override def fileHash = None
        override def lastModifiedTime() = {
          file.map(_.getTimestamp.toInstant).map(FileTime.from).getOrElse({
            throw new FileNotFoundException(s"ftp://$cloudHost/$cloudPath not found")
          })
        }
        override def size() = file.map(_.getSize).getOrElse({
          throw new FileNotFoundException(s"ftp://$cloudHost/$cloudPath not found")
        })
        override def fileKey() = cloudHost + cloudPath
      }
    }
  }

  override def createDirectory(cloudHost: String, cloudPath: String) = {
    runBoolean(FtpOperation(cloudHost, cloudPath, "create directory"), acquireLease(cloudHost))(_.makeDirectory(cloudPath))
  }

  private def findFileSystem(host: String) = getFileSystem(host, fsProvider)

  private def acquireLease[A](host: String): Lease[FTPClient] = findFileSystem(host).leaseClient
}
