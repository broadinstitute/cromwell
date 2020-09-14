package cromwell.filesystems.drs

import cats.effect.IO
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, MarthaResponse}
import common.exception._
import cromwell.core.path.DefaultPathBuilder
import org.apache.commons.lang3.exception.ExceptionUtils
import shapeless.syntax.typeable._


object DrsResolver {
  private val GcsScheme: String = "gs"

  private val GcsProtocolLength: Int = 5 // length of 'gs://'

  def getContainerRelativePath(drsPath: DrsPath): IO[String] = {
    val drsFileSystemProviderOption = drsPath.drsPath.getFileSystem.provider.cast[DrsCloudNioFileSystemProvider]

    val noFileSystemForDrsError = s"Unable to cast file system provider to DrsCloudNioFileSystemProvider for DRS path $drsPath."

    val pathIo = for {
      drsFileSystemProvider <- toIO(drsFileSystemProviderOption, noFileSystemForDrsError)
      marthaResponse <- drsFileSystemProvider.drsPathResolver.resolveDrsThroughMartha(drsPath.pathAsString)
      rootPath = DefaultPathBuilder.get(drsPath.pathWithoutScheme)
      fileName <- getFileName(marthaResponse)
      fullPath = rootPath.resolve(fileName)
      fullPathString = fullPath.pathAsString
    } yield fullPathString


    pathIo
      .handleErrorWith(
        exception =>
          IO.raiseError(
            new RuntimeException(
              s"Error while resolving DRS path: $drsPath. Error: ${ExceptionUtils.getMessage(exception)}"
            )
          )
      )
  }

  /**
    * Return the file name returned from the martha response or get it from the gsUri
    */
  private def getFileName(marthaResponse: MarthaResponse): IO[String] = {
    marthaResponse.fileName match {
      case Some(actualFileName) => IO.pure(actualFileName)
      case None =>
        //Currently, Martha only supports resolving DRS paths to GCS paths
        IO
          .fromEither(marthaResponse.gsUri.toRight(UrlNotFoundException(GcsScheme)))
          .map(_.substring(GcsProtocolLength))
          .map(DefaultPathBuilder.get(_))
          .map(_.name)
    }
  }
}
