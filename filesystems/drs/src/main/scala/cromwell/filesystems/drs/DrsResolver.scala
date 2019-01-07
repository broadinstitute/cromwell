package cromwell.filesystems.drs

import cats.effect.IO
import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, Url}


object DrsResolver {
  private val GcsScheme: String = "gs"

  private def extractPathRelativeToScheme(drsPath: String, urlArray: Array[Url], scheme: String): String = {
    val schemeUrlOption = urlArray.find(u => u.url.startsWith(scheme))

    schemeUrlOption match {
      case Some(schemeUrl) => schemeUrl.url.substring(scheme.length + 3)
      case None => throw UrlNotFoundException(drsPath, scheme)
    }
  }

  def getContainerRelativePath(drsPath: DrsPath): IO[String] = {
    val drsFileSystemProvider = drsPath.drsPath.getFileSystem.provider.asInstanceOf[DrsCloudNioFileSystemProvider]

    for {
      marthaResponse <- drsFileSystemProvider.drsPathResolver.resolveDrsThroughMartha(drsPath.pathAsString)
      //Currently, Martha only supports resolving DRS paths to GCS paths
      relativePath <- IO(extractPathRelativeToScheme(drsPath.pathAsString, marthaResponse.dos.data_object.urls, GcsScheme))
    } yield relativePath
  }
}
