package cloud.nio.impl.ftp

import java.io.IOException
import java.nio.channels.{Channels, ReadableByteChannel, WritableByteChannel}
import java.nio.file.FileAlreadyExistsException

import cats.effect.IO
import cloud.nio.impl.ftp.FtpUtil._
import cloud.nio.impl.ftp.operations._
import cloud.nio.spi.{CloudNioFileList, CloudNioFileProvider, CloudNioRegularFileAttributes}
import cloud.nio.util.TryWithResource._
import io.github.andrebeat.pool.Lease
import org.apache.commons.net.ftp.FTPClient
import org.apache.commons.net.io.Util

import scala.util.{Failure, Success}

class FtpCloudNioFileProvider(fsProvider: FtpCloudNioFileSystemProvider) extends CloudNioFileProvider {
  override def existsPath(cloudHost: String, cloudPath: String): Boolean = withAutoRelease(cloudHost) { client =>
    FtpListFiles(cloudHost, cloudPath, "determine file existence")
      .run(client)
      .map(_.nonEmpty)
  } unsafeRunSync()

  override def existsPaths(cloudHost: String, cloudPathPrefix: String): Boolean = withAutoRelease(cloudHost) { client =>
    existsPathsWithClient(cloudHost, cloudPathPrefix, client)
  } unsafeRunSync()

  private def existsPathsWithClient(cloudHost: String, cloudPathPrefix: String, client: FTPClient): IO[Boolean] = {
    val operation = FtpListDirectories(cloudHost, cloudPathPrefix, "determine directory existence")
    operation
      .run(client)
      .map(_.exists(_.getName == operation.directoryName))
  }

  /**
    * Returns a listing of keys within a bucket starting with prefix. The returned keys should include the prefix. The
    * paths must be absolute, but the key should not begin with a slash.
    */
  override def listObjects(cloudHost: String, cloudPathPrefix: String, markerOption: Option[String]): CloudNioFileList = withAutoRelease(cloudHost) { client =>
    FtpListFiles(cloudHost, cloudPathPrefix, "list objects")
      .run(client)
      .map({ files =>
        val cleanFiles = files.map(_.getName).map(cloudPathPrefix.stripPrefix("/").ensureSlashed + _)
        CloudNioFileList(cleanFiles, markerOption)
      })
  } unsafeRunSync()

  override def copy(sourceCloudHost: String, sourceCloudPath: String, targetCloudHost: String, targetCloudPath: String): Unit = {
    if (sourceCloudHost != targetCloudHost) throw new UnsupportedOperationException(s"Cannot copy files across different ftp servers: Source host: $sourceCloudHost, Target host: $targetCloudHost")

    val fileSystem = findFileSystem(sourceCloudHost)

    val (inputLease, outputLease) = fileSystem.leaseClientsPair

    val streams = () => {
      val is = inputStream(sourceCloudHost, sourceCloudPath, 0, inputLease).unsafeRunSync()
      val os = outputStream(targetCloudHost, targetCloudPath, outputLease).unsafeRunSync()
      new InputOutputStreams(is, os)
    }

    tryWithResource(streams) { ios =>
      Util.copyStream(ios.inputStream, ios.outputStream)
    } match {
      case Success(_) =>
      case Failure(failure) => throw new IOException(s"Failed to copy ftp://$sourceCloudHost/$sourceCloudPath to ftp://$targetCloudHost/$targetCloudPath", failure)
    }
  }

  override def deleteIfExists(cloudHost: String, cloudPath: String): Boolean = withAutoRelease(cloudHost) { client =>
    FtpDeleteFile(cloudHost, cloudPath, "delete").run(client)
  } unsafeRunSync()

  private def inputStream(cloudHost: String, cloudPath: String, offset: Long, lease: Lease[FTPClient]): IO[LeasedInputStream] = {
    FtpInputStream(cloudHost, cloudPath, offset)
      .run(lease.get)
      // Wrap the input stream in a LeasedInputStream so that the lease can be released when the stream is closed
      .map(new LeasedInputStream(cloudHost, cloudPath, _, lease))
  }

  override def read(cloudHost: String, cloudPath: String, offset: Long): ReadableByteChannel = {
    for {
      lease <- acquireLease(cloudHost)
      is <- inputStream(cloudHost, cloudPath, offset, lease)
    } yield Channels.newChannel(is)
  } unsafeRunSync()

  private def outputStream(cloudHost: String, cloudPath: String, lease: Lease[FTPClient]): IO[LeasedOutputStream] = {
    FtpOutputStream(cloudHost, cloudPath)
      .run(lease.get())
      .map(new LeasedOutputStream(cloudHost, cloudPath, _, lease))
  }

  override def write(cloudHost: String, cloudPath: String): WritableByteChannel = {
    for {
      lease <- acquireLease(cloudHost)
      os <- outputStream(cloudHost, cloudPath, lease)
    } yield Channels.newChannel(os)
  } unsafeRunSync()

  override def fileAttributes(cloudHost: String, cloudPath: String): Option[CloudNioRegularFileAttributes] = withAutoRelease(cloudHost) { client =>
    FtpListFiles(cloudHost, cloudPath, "get file attributes")
      .run(client)
      .map(
        _.headOption map { file =>
          new FtpCloudNioRegularFileAttributes(file, cloudHost + cloudPath)
        }
      )
  } unsafeRunSync()

  override def createDirectory(cloudHost: String, cloudPath: String) = withAutoRelease(cloudHost) { client =>
    val operation = FtpCreateDirectory(cloudHost, cloudPath)

    operation.run(client) handleErrorWith {
      /*
       * Sometimes the creation fails with a cryptic error message and the exception generator did not recognize it.
       * In that case, check after the fact if the directory does exist, and if so throw a more appropriate exception 
       */
      case e: FtpIoException =>
        existsPathsWithClient(cloudHost, cloudPath, client) flatMap {
          // If the directory doesn't exist raise the original exception
          case false => IO.raiseError(e)
          // If it does exist, raise a FileAlreadyExistsException
          case _ => IO.raiseError(new FileAlreadyExistsException(operation.fullPath))
        }
      case other => IO.raiseError(other)
    }
  }.void unsafeRunSync()

  private def findFileSystem(host: String): FtpCloudNioFileSystem = fsProvider.newCloudNioFileSystemFromHost(host)

  private def acquireLease[A](host: String): IO[Lease[FTPClient]] = IO { findFileSystem(host).leaseClient }

  private def withAutoRelease[A](cloudHost: String): (FTPClient => IO[A]) => IO[A] = autoRelease[A](acquireLease(cloudHost))
}
