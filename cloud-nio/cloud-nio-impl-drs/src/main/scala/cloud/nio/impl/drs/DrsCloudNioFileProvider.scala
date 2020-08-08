package cloud.nio.impl.drs

import java.nio.channels.{ReadableByteChannel, WritableByteChannel}

import cats.effect.IO
import cloud.nio.spi.{CloudNioFileList, CloudNioFileProvider, CloudNioRegularFileAttributes}
import common.exception._
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.HttpStatus
import org.apache.http.impl.client.HttpClientBuilder


class DrsCloudNioFileProvider(scheme: String,
                              drsPathResolver: EngineDrsPathResolver,
                              httpClientBuilder: HttpClientBuilder,
                              drsReadInterpreter: MarthaResponse => IO[ReadableByteChannel]) extends CloudNioFileProvider {

  private def getDrsPath(cloudHost: String, cloudPath: String): String = s"$scheme://$cloudHost/$cloudPath"

  private def checkIfPathExistsThroughMartha(drsPath: String): IO[Boolean] = {
    /*
     * Unlike other cloud providers where directories are identified with a trailing slash at the end like `gs://bucket/dir/`,
     * DRS has a concept of bundles for directories (not supported yet). Hence for method `checkDirectoryExists` which appends a trailing '/'
     * to see if the current path is a directory, return false
     */
    if (drsPath.endsWith("/")) IO(false)
    else {
      drsPathResolver.rawMarthaResponse(drsPath).use { marthaResponse =>
        val errorMsg = s"Status line was null for martha response $marthaResponse."
        toIO(Option(marthaResponse.getStatusLine), errorMsg)
      }.map(_.getStatusCode == HttpStatus.SC_OK)
    }
  }


  override def existsPath(cloudHost: String, cloudPath: String): Boolean =
    checkIfPathExistsThroughMartha(getDrsPath(cloudHost, cloudPath)).unsafeRunSync()


  override def existsPaths(cloudHost: String, cloudPathPrefix: String): Boolean =
    existsPath(cloudHost, cloudPathPrefix)


  override def listObjects(cloudHost: String, cloudPathPrefix: String, markerOption: Option[String]): CloudNioFileList = {
    val exists = existsPath(cloudHost, cloudPathPrefix)
    val list = if (exists) List(cloudPathPrefix) else Nil
    CloudNioFileList(list, None)
  }


  override def copy(sourceCloudHost: String, sourceCloudPath: String, targetCloudHost: String, targetCloudPath: String): Unit =
    throw new UnsupportedOperationException("DRS currently doesn't support copy.")


  override def deleteIfExists(cloudHost: String, cloudPath: String): Boolean =
    throw new UnsupportedOperationException("DRS currently doesn't support delete.")


  override def read(cloudHost: String, cloudPath: String, offset: Long): ReadableByteChannel = {
    val drsPath = getDrsPath(cloudHost,cloudPath)

    val byteChannelIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath)
      byteChannel <- drsReadInterpreter(marthaResponse)
    } yield byteChannel

    byteChannelIO.handleErrorWith {
        e => IO.raiseError(new RuntimeException(s"Error while reading from DRS path: $drsPath. Error: ${ExceptionUtils.getMessage(e)}"))
    }.unsafeRunSync()
  }


  override def write(cloudHost: String, cloudPath: String): WritableByteChannel =
    throw new UnsupportedOperationException("DRS currently doesn't support write.")


  override def fileAttributes(cloudHost: String, cloudPath: String): Option[CloudNioRegularFileAttributes] =
    Option(new DrsCloudNioRegularFileAttributes(getDrsPath(cloudHost,cloudPath), drsPathResolver))
}



case class GcsFilePath(bucket: String, file: String)

