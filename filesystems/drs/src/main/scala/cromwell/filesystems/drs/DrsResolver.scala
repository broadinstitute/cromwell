package cromwell.filesystems.drs

import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, MarthaResponse, CloudUrl}


object DrsResolver {
  private val GcsScheme: String = "gs"

  private def extractPathRelativeToScheme(drsPath: String, urlArray: Array[CloudUrl], scheme: String): String = {
    val schemeUrlOption = urlArray.find(u => u.url.startsWith(scheme))

    schemeUrlOption match {
      case Some(schemeUrl) => schemeUrl.url.substring(scheme.length + 3)
      case None => throw UrlNotFoundException(drsPath, scheme)
    }
  }

  def getContainerRelativePath(drsPath: DrsPath): String = {
    val drsFileSystemProvider = drsPath.drsPath.getFileSystem.provider.asInstanceOf[DrsCloudNioFileSystemProvider]
    val marthaResponseObj: MarthaResponse = drsFileSystemProvider.drsPathResolver.resolveDrsThroughMartha(drsPath.pathAsString)

    //Currently, Martha only supports resolving DRS paths to GCS paths
    extractPathRelativeToScheme(drsPath.pathAsString, marthaResponseObj.dos.data_object.urls, GcsScheme)
  }
}


case class UrlNotFoundException(drsPath: String, scheme: String) extends Exception(s"DRS was not able to find a $scheme url associated with $drsPath.")
