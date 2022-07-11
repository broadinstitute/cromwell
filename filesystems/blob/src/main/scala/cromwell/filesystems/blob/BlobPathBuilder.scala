package cromwell.filesystems.blob

import cromwell.core.path.{NioPath, PathBuilder}
import scala.util.Try
import java.nio.file.FileSystems
import java.net.URI
import scala.language.postfixOps
import com.google.common.net.UrlEscapers
import com.azure.storage.blob.nio.AzureFileSystem
import scala.jdk.CollectionConverters._
import cromwell.core.path.Path
import scala.util.Failure
import com.azure.core.credential.AzureSasCredential
import cromwell.filesystems.blob.BlobPathBuilder._

object BlobPathBuilder {

  sealed trait BlobPathValidation
  case class ValidBlobPath(path: String) extends BlobPathValidation
  case class UnparsableBlobPath(errorMessage: Throwable) extends BlobPathValidation

}

class BlobPathBuilder(credential: AzureSasCredential, container: String, storageAccount: String) extends PathBuilder {

  val fileSystemConfig: Map[String, Object] = Map((AzureFileSystem.AZURE_STORAGE_SAS_TOKEN_CREDENTIAL, credential),
                                                  (AzureFileSystem.AZURE_STORAGE_FILE_STORES, container))

  def validateBlobPath(string: String): BlobPathValidation = {
    Try {
      val uri = URI.create(UrlEscapers.urlFragmentEscaper().escape(string))
      ValidBlobPath(uri.getPath.replaceFirst("/" + container, ""))
    } recover{ case t => UnparsableBlobPath(t) } get
  }

  def endpointString = "https://" + storageAccount + ".blob.core.windows.net"

  def build(string: String): Try[BlobPath] = {
    validateBlobPath(string) match {
      case ValidBlobPath(path) =>
        Try {
          val fileSystem = FileSystems.newFileSystem(new URI("azb://?endpoint=" + endpointString), fileSystemConfig.asJava)
          val blobStoragePath = fileSystem.getPath(path)
          BlobPath(blobStoragePath, endpointString, container)
        }
      case UnparsableBlobPath(errorMessage: Throwable) => Failure(errorMessage)
    }
  }

  override def name: String = "Azure Blob Storage"
}

// Add args for container, storage account name
case class BlobPath private[blob](nioPath: NioPath, endpoint: String, container: String) extends Path {
  override protected def newPath(nioPath: NioPath): Path = BlobPath(nioPath, endpoint, container)

  override def pathAsString: String = endpoint + "/" + container + "/" + nioPath.toString()

  override def pathWithoutScheme: String = "/" + container + "/" + nioPath.toString()
}