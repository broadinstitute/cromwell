package cromwell.filesystems.drs

import cats.data.NonEmptyList
import cats.effect.IO
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, DrsPathResolver, MarthaField}
import common.exception._
import cromwell.core.path.{DefaultPathBuilder, Path}
import org.apache.commons.lang3.exception.ExceptionUtils
import shapeless.syntax.typeable._


object DrsResolver {
  private val GcsScheme: String = "gs"

  private val GcsProtocolLength: Int = 5 // length of 'gs://'

  private def resolveError[A](pathAsString: String)(throwable: Throwable): IO[A] = {
    IO.raiseError(
      new RuntimeException(
        s"Error while resolving DRS path: $pathAsString. Error: ${ExceptionUtils.getMessage(throwable)}"
      )
    )
  }

  private def getDrsPathResolver(drsPath: DrsPath): IO[DrsPathResolver] = {
    val drsFileSystemProviderOption = drsPath.drsPath.getFileSystem.provider.cast[DrsCloudNioFileSystemProvider]

    val noFileSystemForDrsError = s"Unable to cast file system provider to DrsCloudNioFileSystemProvider for DRS path $drsPath."

    for {
      drsFileSystemProvider <- toIO(drsFileSystemProviderOption, noFileSystemForDrsError)
    } yield drsFileSystemProvider.drsPathResolver
  }

  case class MarthaLocalizationData(gsUri: Option[String],
                                    fileName: Option[String],
                                    bondProvider: Option[String],
                                    localizationPath: Option[String])

  private def getMarthaLocalizationData(pathAsString: String,
                                        drsPathResolver: DrsPathResolver): IO[MarthaLocalizationData] = {
    val fields = NonEmptyList.of(MarthaField.GsUri, MarthaField.FileName, MarthaField.BondProvider, MarthaField.LocalizationPath)

    drsPathResolver.resolveDrsThroughMartha(pathAsString, fields) map { r =>
      MarthaLocalizationData(r.gsUri, r.fileName, r.bondProvider, r.localizationPath)
    }
  }

  /** Returns the `gsUri` if it ends in the `fileName` and the `bondProvider` is empty. */
  private def getSimpleGsUri(localizationData: MarthaLocalizationData): Option[String] = {
    localizationData match {
      // `gsUri` not defined so no gsUri can be returned.
      case MarthaLocalizationData(None, _, _, _) => None
      // `bondProvider` defined, cannot "preresolve" to GCS.
      case MarthaLocalizationData(_, _, Some(_), _) => None
      // `localizationPath` defined which takes precedence over `fileName`. Don't attempt preresolve for this case.
      case MarthaLocalizationData(_, _, _ , Some(_)) => None
      case MarthaLocalizationData(Some(gsUri), Some(fileName), _, _) if !gsUri.endsWith(s"/$fileName") => None
      case MarthaLocalizationData(Some(gsUri), _, _, _) => Option(gsUri)
    }
  }

  /** Returns the `gsUri` if it ends in the `fileName` and the `bondProvider` is empty. */
  def getSimpleGsUri(pathAsString: String,
                     drsPathResolver: DrsPathResolver): IO[Option[String]] = {

    val gsUriIO = getMarthaLocalizationData(pathAsString, drsPathResolver) map getSimpleGsUri

    gsUriIO.handleErrorWith(resolveError(pathAsString))
  }

  /** Returns the `gsUri` if it ends in the `fileName` and the `bondProvider` is empty. */
  def getSimpleGsUri(drsPath: DrsPath): IO[Option[String]] = {
    for {
      drsPathResolver <- getDrsPathResolver(drsPath)
      gsUri <- getSimpleGsUri(drsPath.pathAsString, drsPathResolver)
    } yield gsUri
  }

  def getContainerRelativePath(drsPath: DrsPath): IO[String] = {

    def buildFullPath(fileName: String, rootPath: Path): IO[Path] = {
      if (fileName.startsWith("/"))
        // If `fileName` has a leading slash do not treat it as being relative to the DRS `rootPath` as this is the way
        // localization paths are specified by TDR [BT-418]. Do strip the leading slash as the resulting path will still
        // need to be made relative to the container root.
        IO.fromTry(DefaultPathBuilder.build(fileName.substring(1)))
      else {
        // Regular non-leading-slash `fileName` should be made relative to the `rootPath`.
        IO(rootPath.resolve(fileName))
      }
    }

    val pathIO = for {
      drsPathResolver <- getDrsPathResolver(drsPath)
      localizationData <- getMarthaLocalizationData(drsPath.pathAsString, drsPathResolver)
      /*
      In the DOS/DRS spec file names are safe for file systems but not necessarily the DRS URIs.
      Reuse the regex defined for ContentsObject.name, plus add "/" for directory separators.
      https://ga4gh.github.io/data-repository-service-schemas/preview/release/drs-1.0.0/docs/#_contentsobject
       */
      rootPath = DefaultPathBuilder.get(drsPath.pathWithoutScheme.replaceAll("[^/A-Za-z0-9._-]", "_"))
      fileName <- getFileName(localizationData)
      fullPath <- buildFullPath(fileName, rootPath)
    } yield fullPath.pathAsString

    pathIO.handleErrorWith(resolveError(drsPath.pathAsString))
  }

  /**
    * Return the file name returned from the martha response or get it from the gsUri
    */
  private def getFileName(localizationData: MarthaLocalizationData): IO[String] = {
    localizationData match {
      case MarthaLocalizationData(_, _, _, Some(localizationPath)) => IO.pure(localizationPath)
      case MarthaLocalizationData(_, Some(fileName), _, _) => IO.pure(fileName)
      case _ =>
        // If this logic is forced to fall back on the GCS path there better be a GCS path to fall back on.
        IO
          .fromEither(localizationData.gsUri.toRight(UrlNotFoundException(GcsScheme)))
          .map(_.substring(GcsProtocolLength))
          .map(DefaultPathBuilder.get(_))
          .map(_.name)
    }
  }
}
