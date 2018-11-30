package cromwell.filesystems.drs

import cloud.nio.impl.drs.{DrsPathResolver, MarthaResponse, Url}
import com.typesafe.scalalogging.StrictLogging


class DrsResolver extends StrictLogging {
  private val GcsScheme: String = "gs"

  private def extractPathRelativeToScheme(drsPath: String, urlArray: Array[Url], scheme: String): String = {
    val schemeUrlOption = urlArray.find(u => u.url.startsWith(scheme))

    schemeUrlOption match {
      case Some(schemeUrl) => schemeUrl.url.substring(scheme.length + 3)
      case None => throw GcsUrlNotFoundException(drsPath, scheme)
    }
  }

  def getContainerRelativePath(drsPath: DrsPath): String = {
    val drsFileSystemGlobalConfig = drsPath.drsPath.filesystem.provider.config
    val drsPathResolver: DrsPathResolver = DrsPathResolver(drsFileSystemGlobalConfig)

    val marthaResponseObj: MarthaResponse = drsPathResolver.resolveDrsThroughMartha(drsPath.pathAsString)

    //Currently, Martha only supports resolving DRS paths to GCS paths
    extractPathRelativeToScheme(drsPath.pathAsString, marthaResponseObj.dos.data_object.urls, GcsScheme)
  }
}


case class GcsUrlNotFoundException(drsPath: String, scheme: String) extends Exception(s"DRS was not able to find a $scheme url associated with $drsPath.")
