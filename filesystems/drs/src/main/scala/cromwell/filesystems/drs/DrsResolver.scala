package cromwell.filesystems.drs

import cats.data.NonEmptyList
import cats.effect.IO
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, DrsPathResolver, DrsResolverField}
import common.exception._
import cromwell.core.path.{DefaultPathBuilder, Path}
import org.apache.commons.lang3.exception.ExceptionUtils
import shapeless.syntax.typeable._

object DrsResolver {
  private val GcsScheme: String = "gs"

  private val GcsProtocolLength: Int = 5 // length of 'gs://'

  private val DrsLocalizationPathsContainer = "drs_localization_paths"

  private def resolveError[A](pathAsString: String)(throwable: Throwable): IO[A] =
    IO.raiseError(
      new RuntimeException(
        s"Error while resolving DRS path: $pathAsString. Error: ${ExceptionUtils.getMessage(throwable)}"
      )
    )

  private def getDrsPathResolver(drsPath: DrsPath): IO[DrsPathResolver] = {
    val drsFileSystemProviderOption = drsPath.drsPath.getFileSystem.provider.cast[DrsCloudNioFileSystemProvider]

    val noFileSystemForDrsError =
      s"Unable to cast file system provider to DrsCloudNioFileSystemProvider for DRS path $drsPath."

    for {
      drsFileSystemProvider <- toIO(drsFileSystemProviderOption, noFileSystemForDrsError)
    } yield drsFileSystemProvider.drsPathResolver
  }

  case class DrsResolverLocalizationData(gsUri: Option[String],
                                         fileName: Option[String],
                                         bondProvider: Option[String],
                                         localizationPath: Option[String]
  )

  private def getDrsResolverLocalizationData(pathAsString: String,
                                             drsPathResolver: DrsPathResolver
  ): IO[DrsResolverLocalizationData] = {
    val fields = NonEmptyList.of(DrsResolverField.GsUri,
                                 DrsResolverField.FileName,
                                 DrsResolverField.BondProvider,
                                 DrsResolverField.LocalizationPath
    )

    drsPathResolver.resolveDrs(pathAsString, fields) map { r =>
      DrsResolverLocalizationData(r.gsUri, r.fileName, r.bondProvider, r.localizationPath)
    }
  }

  /** Returns the `gsUri` if it ends in the `fileName` and the `bondProvider` is empty. */
  private def getSimpleGsUri(localizationData: DrsResolverLocalizationData): Option[String] =
    localizationData match {
      // `gsUri` not defined so no gsUri can be returned.
      case DrsResolverLocalizationData(None, _, _, _) => None
      // `bondProvider` defined, cannot "preresolve" to GCS.
      case DrsResolverLocalizationData(_, _, Some(_), _) => None
      // Do not return the simple GS URI if the `fileName` from metadata is mismatched to the filename in the `gsUri`.
      case DrsResolverLocalizationData(Some(gsUri), Some(fileName), _, _) if !gsUri.endsWith(s"/$fileName") => None
      // Barring any of the situations above return the `gsUri`.
      case DrsResolverLocalizationData(Some(gsUri), _, _, _) => Option(gsUri)
    }

  /** Returns the `gsUri` if it ends in the `fileName` and the `bondProvider` is empty. */
  def getSimpleGsUri(pathAsString: String, drsPathResolver: DrsPathResolver): IO[Option[String]] = {

    val gsUriIO = getDrsResolverLocalizationData(pathAsString, drsPathResolver) map getSimpleGsUri

    gsUriIO.handleErrorWith(resolveError(pathAsString))
  }

  /** Returns the `gsUri` if it ends in the `fileName` and the `bondProvider` is empty. */
  def getSimpleGsUri(drsPath: DrsPath): IO[Option[String]] =
    for {
      drsPathResolver <- getDrsPathResolver(drsPath)
      gsUri <- getSimpleGsUri(drsPath.pathAsString, drsPathResolver)
    } yield gsUri

  def getContainerRelativePath(drsPath: DrsPath): IO[String] = {
    val pathIO = for {
      drsPathResolver <- getDrsPathResolver(drsPath)
      localizationData <- getDrsResolverLocalizationData(drsPath.pathAsString, drsPathResolver)
      containerRelativePath <- buildContainerRelativePath(localizationData, drsPath)
    } yield containerRelativePath.pathAsString

    pathIO.handleErrorWith(resolveError(drsPath.pathAsString))
  }

  // Return the container relative path built from the DRS Resolver-specified localization path, file name, or gs URI.
  private def buildContainerRelativePath(localizationData: DrsResolverLocalizationData, drsPath: Path): IO[Path] = {
    // Return a relative path constructed from the DRS path minus the leading scheme.
    // In the DOS/DRS spec file names are safe for file systems but not necessarily the DRS URIs.
    // Reuse the regex defined for ContentsObject.name, plus add "/" for directory separators.
    // https://ga4gh.github.io/data-repository-service-schemas/preview/release/drs-1.0.0/docs/#_contentsobject
    def drsPathRelativePath: Path =
      DefaultPathBuilder.get(drsPath.pathWithoutScheme.replaceAll("[^/A-Za-z0-9._-]", "_"))

    localizationData match {
      case DrsResolverLocalizationData(_, _, _, Some(localizationPath)) =>
        // TDR may return an explicit localization path and if so we should not use the `drsPathRelativePath`.
        // We want to end up with something like /cromwell_root/drs_localization_paths/tdr/specified/path/foo.bam.
        // Calling code will add the `/cromwell_root/`, so strip any leading slashes to make this a relative path:
        val relativeLocalizationPath = if (localizationPath.startsWith("/")) localizationPath.tail else localizationPath
        IO.fromTry(DefaultPathBuilder.build(DrsLocalizationPathsContainer).map(_.resolve(relativeLocalizationPath)))
      case DrsResolverLocalizationData(_, Some(fileName), _, _) =>
        // Paths specified by filename only are made relative to `drsPathRelativePath`.
        IO(drsPathRelativePath.resolve(fileName))
      case _ =>
        // If this logic is forced to fall back on the GCS path there better be a GCS path to fall back on.
        IO
          .fromEither(localizationData.gsUri.toRight(UrlNotFoundException(GcsScheme)))
          .map(_.substring(GcsProtocolLength))
          .map(DefaultPathBuilder.get(_))
          .map(path => drsPathRelativePath.resolve(path.name))
    }
  }
}
