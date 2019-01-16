package cromwell.filesystems.drs

import cats.effect.IO
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, Url}
import common.exception._
import org.apache.commons.lang3.exception.ExceptionUtils
import shapeless.syntax.typeable._


object DrsResolver {
  private val GcsScheme: String = "gs"

  private def urlProtocolLength(scheme: String): Int = scheme.length + 3 //3: length of '://'


  def extractUrlForScheme(urlArray: Array[Url], scheme: String): Either[UrlNotFoundException, String] = {
    val schemeUrlOption = urlArray.find(_.url.startsWith(scheme))

    schemeUrlOption match {
      case Some(schemeUrl) => Right(schemeUrl.url)
      case None => Left(UrlNotFoundException(scheme))
    }
  }

  def getContainerRelativePath(drsPath: DrsPath): IO[String] = {
    val drsFileSystemProviderOption = drsPath.drsPath.getFileSystem.provider.cast[DrsCloudNioFileSystemProvider]

    val noFileSystemForDrsError = s"Unable to cast file system provider to DrsCloudNioFileSystemProvider for DRS path $drsPath."

    for {
      drsFileSystemProvider <- toIO(drsFileSystemProviderOption, noFileSystemForDrsError)
      marthaResponse <- drsFileSystemProvider.drsPathResolver.resolveDrsThroughMartha(drsPath.pathAsString)
      //Currently, Martha only supports resolving DRS paths to GCS paths
      relativePath <- IO.fromEither(extractUrlForScheme(marthaResponse.dos.data_object.urls, GcsScheme))
        .map(_.substring(urlProtocolLength(GcsScheme)))
        .handleErrorWith(e => IO.raiseError(new RuntimeException(s"Error while resolving DRS path: $drsPath. Error: ${ExceptionUtils.getMessage(e)}")))
    } yield relativePath
  }
}
