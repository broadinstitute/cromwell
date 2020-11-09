package cromwell.filesystems.drs

import cats.data.NonEmptyList
import cats.effect.IO
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, MarthaField}
import common.exception._
import cromwell.core.path.DefaultPathBuilder
import org.apache.commons.lang3.exception.ExceptionUtils
import shapeless.syntax.typeable._


object DrsResolver {
  private val GcsScheme: String = "gs"

  private val GcsProtocolLength: Int = 5 // length of 'gs://'

  private def resolveError[A](drsPath: DrsPath)(throwable: Throwable): IO[A] = {
    IO.raiseError(
      new RuntimeException(
        s"Error while resolving DRS path: $drsPath. Error: ${ExceptionUtils.getMessage(throwable)}"
      )
    )
  }

  private def getFileNameAndGsUri(drsPath: DrsPath): IO[(Option[String], Option[String])] = {
    val drsFileSystemProviderOption = drsPath.drsPath.getFileSystem.provider.cast[DrsCloudNioFileSystemProvider]

    val noFileSystemForDrsError = s"Unable to cast file system provider to DrsCloudNioFileSystemProvider for DRS path $drsPath."
    val fields = NonEmptyList.of(MarthaField.FileName, MarthaField.GsUri)

    for {
      drsFileSystemProvider <- toIO(drsFileSystemProviderOption, noFileSystemForDrsError)
      marthaResponse <- drsFileSystemProvider.drsPathResolver.resolveDrsThroughMartha(drsPath.pathAsString, fields)
    } yield (marthaResponse.fileName, marthaResponse.gsUri)
  }

  def getGsUri(drsPath: DrsPath): IO[String] = {
    val gsUriIO = for {
      fileNameAndGsUri <- getFileNameAndGsUri(drsPath)
      (_, gsUriOption) = fileNameAndGsUri
      gsUri <- IO.fromEither(gsUriOption.toRight(UrlNotFoundException(GcsScheme)))
    } yield gsUri

    gsUriIO.handleErrorWith(resolveError(drsPath))
  }

  def getContainerRelativePath(drsPath: DrsPath): IO[String] = {
    val pathIO = for {
      fileNameAndGsUri <- getFileNameAndGsUri(drsPath)
      (fileNameOption, gsUriOption) = fileNameAndGsUri
      /*
      In the DOS/DRS spec file names are safe for file systems but not necessarily the DRS URIs.
      Reuse the regex defined for ContentsObject.name, plus add "/" for directory separators.
      https://ga4gh.github.io/data-repository-service-schemas/preview/release/drs-1.0.0/docs/#_contentsobject
       */
      rootPath = DefaultPathBuilder.get(drsPath.pathWithoutScheme.replaceAll("[^/A-Za-z0-9._-]", "_"))
      fileName <- getFileName(fileNameOption, gsUriOption)
      fullPath = rootPath.resolve(fileName)
      fullPathString = fullPath.pathAsString
    } yield fullPathString

    pathIO.handleErrorWith(resolveError(drsPath))
  }

  /**
    * Return the file name returned from the martha response or get it from the gsUri
    */
  private def getFileName(fileName: Option[String], gsUri: Option[String]): IO[String] = {
    fileName match {
      case Some(actualFileName) => IO.pure(actualFileName)
      case None =>
        //Currently, Martha only supports resolving DRS paths to GCS paths
        IO
          .fromEither(gsUri.toRight(UrlNotFoundException(GcsScheme)))
          .map(_.substring(GcsProtocolLength))
          .map(DefaultPathBuilder.get(_))
          .map(_.name)
    }
  }
}
