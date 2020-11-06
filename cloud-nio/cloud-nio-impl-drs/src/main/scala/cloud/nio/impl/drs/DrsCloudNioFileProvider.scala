package cloud.nio.impl.drs

import java.nio.channels.{ReadableByteChannel, WritableByteChannel}

import cats.data.NonEmptyList
import cats.effect.IO
import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.impl.drs.DrsCloudNioRegularFileAttributes._
import cloud.nio.spi.{CloudNioFileList, CloudNioFileProvider, CloudNioRegularFileAttributes}
import common.exception._
import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.http.HttpStatus


class DrsCloudNioFileProvider(drsPathResolver: EngineDrsPathResolver,
                              drsReadInterpreter: DrsReadInterpreter) extends CloudNioFileProvider {

  private def checkIfPathExistsThroughMartha(drsPath: String): IO[Boolean] = {
    /*
     * Unlike other cloud providers where directories are identified with a trailing slash at the end like `gs://bucket/dir/`,
     * DRS has a concept of bundles for directories (not supported yet). Hence for method `checkDirectoryExists` which appends a trailing '/'
     * to see if the current path is a directory, return false
     */
    if (drsPath.endsWith("/")) IO(false)
    else {
      drsPathResolver.rawMarthaResponse(drsPath, NonEmptyList.one(MarthaField.GsUri)).use { marthaResponse =>
        val errorMsg = s"Status line was null for martha response $marthaResponse."
        toIO(Option(marthaResponse.getStatusLine), errorMsg)
      }.map(_.getStatusCode == HttpStatus.SC_OK)
    }
  }

  override def existsPath(drsPath: String, unused: String): Boolean =
    checkIfPathExistsThroughMartha(drsPath).unsafeRunSync()

  override def existsPaths(cloudHost: String, cloudPathPrefix: String): Boolean =
    existsPath(cloudHost, cloudPathPrefix)

  override def listObjects(drsPath: String, unused: String, markerOption: Option[String]): CloudNioFileList = {
    throw new UnsupportedOperationException("DRS currently doesn't support list.")
  }

  override def copy(sourceCloudHost: String, sourceCloudPath: String, targetCloudHost: String, targetCloudPath: String): Unit =
    throw new UnsupportedOperationException("DRS currently doesn't support copy.")

  override def deleteIfExists(cloudHost: String, cloudPath: String): Boolean =
    throw new UnsupportedOperationException("DRS currently doesn't support delete.")

  override def read(drsPath: String, unused: String, offset: Long): ReadableByteChannel = {
    val fields = NonEmptyList.of(MarthaField.GsUri, MarthaField.GoogleServiceAccount)

    val byteChannelIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath, fields)
      byteChannel <- drsReadInterpreter(marthaResponse.gsUri, marthaResponse.googleServiceAccount)
    } yield byteChannel

    byteChannelIO.handleErrorWith {
        e => IO.raiseError(new RuntimeException(s"Error while reading from DRS path: $drsPath. Error: ${ExceptionUtils.getMessage(e)}"))
    }.unsafeRunSync()
  }

  override def write(cloudHost: String, cloudPath: String): WritableByteChannel =
    throw new UnsupportedOperationException("DRS currently doesn't support write.")

  override def fileAttributes(drsPath: String, unused: String): Option[CloudNioRegularFileAttributes] = {
    val fields = NonEmptyList.of(MarthaField.Size, MarthaField.TimeCreated, MarthaField.TimeUpdated, MarthaField.Hashes)

    val fileAttributesIO = for {
      marthaResponse <- drsPathResolver.resolveDrsThroughMartha(drsPath, fields)
      sizeOption = marthaResponse.size
      hashoption = getPreferredHash(marthaResponse.hashes)
      timeCreatedOption <- convertToFileTime(drsPath, MarthaField.TimeCreated, marthaResponse.timeCreated)
      timeUpdatedOption <- convertToFileTime(drsPath, MarthaField.TimeUpdated, marthaResponse.timeUpdated)
    } yield new DrsCloudNioRegularFileAttributes(drsPath, sizeOption, hashoption, timeCreatedOption, timeUpdatedOption)

    Option(fileAttributesIO.unsafeRunSync())
  }
}

object DrsCloudNioFileProvider {
  type DrsReadInterpreter = (Option[String], Option[SADataObject]) => IO[ReadableByteChannel]
}
