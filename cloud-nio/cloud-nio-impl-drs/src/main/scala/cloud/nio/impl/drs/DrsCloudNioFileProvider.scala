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

  private def checkIfPathExistsThroughDrsResolver(drsPath: String): IO[Boolean] = {
    /*
     * Unlike other cloud providers where directories are identified with a trailing slash at the end like `gs://bucket/dir/`,
     * DRS has a concept of bundles for directories (not supported yet). Hence for method `checkDirectoryExists` which appends a trailing '/'
     * to see if the current path is a directory, return false
     */
    if (drsPath.endsWith("/")) IO(false)
    else {
      drsPathResolver.rawDrsResolverResponse(drsPath, NonEmptyList.one(DrsResolverField.GsUri)).use { drsResolverResponse =>
        val errorMsg = s"Status line was null for DRS Resolver response $drsResolverResponse."
        toIO(Option(drsResolverResponse.getStatusLine), errorMsg)
      }.map(_.getStatusCode == HttpStatus.SC_OK)
    }
  }

  override def existsPath(drsPath: String, unused: String): Boolean =
    checkIfPathExistsThroughDrsResolver(drsPath).unsafeRunSync()

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
    val fields = NonEmptyList.of(DrsResolverField.GsUri, DrsResolverField.GoogleServiceAccount, DrsResolverField.AccessUrl)

    val byteChannelIO = for {
      drsResolverResponse <- drsPathResolver.resolveDrs(drsPath, fields)
      byteChannel <- drsReadInterpreter(drsPathResolver, drsResolverResponse)
    } yield byteChannel

    byteChannelIO.handleErrorWith {
        e => IO.raiseError(new RuntimeException(s"Error while reading from DRS path: $drsPath. Error: ${ExceptionUtils.getMessage(e)}"))
    }.unsafeRunSync()
  }

  override def write(cloudHost: String, cloudPath: String): WritableByteChannel =
    throw new UnsupportedOperationException("DRS currently doesn't support write.")

  override def fileAttributes(drsPath: String, unused: String): Option[CloudNioRegularFileAttributes] = {
    val fields = NonEmptyList.of(DrsResolverField.Size, DrsResolverField.TimeCreated, DrsResolverField.TimeUpdated, DrsResolverField.Hashes)

    val fileAttributesIO = for {
      drsResolverResponse <- drsPathResolver.resolveDrs(drsPath, fields)
      sizeOption = drsResolverResponse.size
      hashOption = getPreferredHash(drsResolverResponse.hashes)
      timeCreatedOption <- convertToFileTime(drsPath, DrsResolverField.TimeCreated, drsResolverResponse.timeCreated)
      timeUpdatedOption <- convertToFileTime(drsPath, DrsResolverField.TimeUpdated, drsResolverResponse.timeUpdated)
    } yield new DrsCloudNioRegularFileAttributes(drsPath, sizeOption, hashOption, timeCreatedOption, timeUpdatedOption)

    Option(fileAttributesIO.unsafeRunSync())
  }
}

object DrsCloudNioFileProvider {
  type DrsReadInterpreter = (DrsPathResolver, DrsResolverResponse) => IO[ReadableByteChannel]
}
