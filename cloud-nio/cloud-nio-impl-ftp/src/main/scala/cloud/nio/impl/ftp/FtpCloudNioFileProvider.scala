package cloud.nio.impl.ftp

import java.io.{IOException, InputStream, OutputStream}
import java.nio.channels.Channels

import cloud.nio.impl.ftp.FtpCloudNioFileProvider._
import cloud.nio.impl.ftp.FtpUtil._
import cloud.nio.spi.{CloudNioFileList, CloudNioFileProvider, CloudNioFileSystem, CloudNioRegularFileAttributes}
import cloud.nio.util.TryWithResource._
import io.github.andrebeat.pool.Lease
import org.apache.commons.net.ftp.FTPClient
import org.apache.commons.net.io.Util

import scala.util.{Failure, Success}

object FtpCloudNioFileProvider {
  implicit class EnhancedCloudPath(val cloudPath: String) extends AnyVal {
    def ensureSlashedPrefix = if (cloudPath.startsWith("/")) cloudPath else s"/$cloudPath"
    def ensureSlashed = if (cloudPath.endsWith("/")) cloudPath else s"$cloudPath/"
  }
}

class FtpCloudNioFileProvider(fsProvider: FtpCloudNioFileSystemProvider) extends CloudNioFileProvider {
  override def existsPath(cloudHost: String, cloudPath: String) = {
    run(FtpOperation(cloudHost, cloudPath, "determine file existence"), acquireLease(cloudHost))(_.listFiles(cloudPath.ensureSlashedPrefix)).nonEmpty
  }

  override def existsPaths(cloudHost: String, cloudPathPrefix: String) = {
    // We need to list the directories in the parent and see if any matches the name, hence the string manipulations
    val parts = cloudPathPrefix.ensureSlashedPrefix.stripSuffix(CloudNioFileSystem.Separator).split(CloudNioFileSystem.Separator)
    val parent = parts.init.mkString(CloudNioFileSystem.Separator)
    val directoryName = parts.last

    run(FtpOperation(cloudHost, cloudPathPrefix, "determine directory existence"), acquireLease(cloudHost))(_.listDirectories(parent))
      .exists(_.getName == directoryName)
  }

  /**
    * Returns a listing of keys within a bucket starting with prefix. The returned keys should include the prefix. The
    * paths must be absolute, but the key should not begin with a slash.
    */
  override def listObjects(cloudHost: String, cloudPathPrefix: String, markerOption: Option[String]): CloudNioFileList = {
    val files = run(FtpOperation(cloudHost, cloudPathPrefix, "list objects"), acquireLease(cloudHost)) { client =>
      client.listFiles(cloudPathPrefix.ensureSlashedPrefix)
    }

    CloudNioFileList(files.map(_.getName).map(cloudPathPrefix.stripPrefix("/").ensureSlashed + _), markerOption)
  }

  override def copy(sourceCloudHost: String, sourceCloudPath: String, targetCloudHost: String, targetCloudPath: String) = {
    if (sourceCloudHost != targetCloudHost) throw new UnsupportedOperationException(s"Cannot copy files across different ftp servers: Source host: $sourceCloudHost, Target host: $targetCloudHost")

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
    runBoolean(FtpOperation(cloudHost, cloudPath, "delete"), acquireLease(cloudHost), failOnFalse = false)(_.deleteFile(cloudPath.ensureSlashedPrefix))
  }

  private def inputStream(cloudHost: String, cloudPath: String, offset: Long, lease: Lease[FTPClient]) = {
    val operation = FtpOperation(cloudHost, cloudPath, "read")

    val stream: InputStream = run(operation, lease, autoRelease = false) { client =>
      client.setRestartOffset(offset)
      client.retrieveFileStream(cloudPath.ensureSlashedPrefix)
    }

    // Wrap the input stream in a LeasedInputStream so that the lease can be release when the stream is closed
    new LeasedInputStream(cloudHost, cloudPath, stream, lease)
  }

  override def read(cloudHost: String, cloudPath: String, offset: Long) = {
    Channels.newChannel(inputStream(cloudHost, cloudPath, offset, acquireLease(cloudHost)))
  }

  private def outputStream(cloudHost: String, cloudPath: String, lease: Lease[FTPClient]) = {
    val operation = FtpOperation(cloudHost, cloudPath, "write")

    val stream: OutputStream = run(operation, lease, autoRelease = false)(_.storeFileStream(cloudPath.ensureSlashedPrefix))

    // Wrap the input stream in a LeasedOutputStream so that the lease can be release when the stream is closed
    new LeasedOutputStream(cloudHost: String, cloudPath: String, stream, lease)
  }

  override def write(cloudHost: String, cloudPath: String) = {
    Channels.newChannel(outputStream(cloudHost: String, cloudPath: String, acquireLease(cloudHost)))
  }

  override def fileAttributes(cloudHost: String, cloudPath: String): Option[CloudNioRegularFileAttributes] = {
    run(FtpOperation(cloudHost, cloudPath, "get file attributes"), acquireLease(cloudHost))(_.listFiles(cloudPath.ensureSlashedPrefix)).headOption map { file =>
      new FtpCloudNioRegularFileAttributes(file, cloudHost + cloudPath)
    }
  }

  override def createDirectory(cloudHost: String, cloudPath: String) = {
    runBoolean(FtpOperation(cloudHost, cloudPath, "create directory"), acquireLease(cloudHost))(_.makeDirectory(cloudPath))
    ()
  }

  private def findFileSystem(host: String) = fsProvider.newCloudNioFileSystemFromHost(host)

  private def acquireLease[A](host: String): Lease[FTPClient] = findFileSystem(host).leaseClient
}
