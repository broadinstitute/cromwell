package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import com.azure.storage.blob.nio.AzureFileSystem
import com.google.common.net.UrlEscapers
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.blob.BlobPathBuilder._

import java.net.{MalformedURLException, URI}
import java.nio.file.{FileSystem, FileSystemNotFoundException, FileSystems}
import scala.jdk.CollectionConverters._
import scala.language.postfixOps
import scala.util.{Failure, Try}

object BlobPathBuilder {

  sealed trait BlobPathValidation
  case class ValidBlobPath(path: String) extends BlobPathValidation
  case class UnparsableBlobPath(errorMessage: Throwable) extends BlobPathValidation

  def invalidBlobPathMessage(container: String, endpoint: String) = s"Malformed Blob URL for this builder. Expecting a URL for a container $container and endpoint $endpoint"
  def parseURI(string: String) = URI.create(UrlEscapers.urlFragmentEscaper().escape(string))
  def parseStorageAccount(uri: URI) = uri.getHost().split("\\.").filter(!_.isEmpty()).headOption

  /**
    * Validates a that a path from a string is a valid BlobPath of the format:
    * {endpoint}/{containerName}/{pathToFile}
    *
    * with an endpoint for a particular storage account typically given by:
    * https://{storageAccountName}.blob.core.windows.net/
    *
    * For example, a path string we might expect to receive might look like:
    * https://appexternalstorage.blob.core.windows.net/inputs/test/testFile.wdl
    *
    * In this example
    * storageAccountName -> appexternalstorage
    * endpoint -> https://{storageAccountName}.blob.core.windows.net/
    * container -> inputs
    * pathToFile -> test/testFile.wdl
    *
    * If the configured container and storage account do not match, the string is considered unparsable
    */
  def validateBlobPath(string: String, container: String, endpoint: String): BlobPathValidation = {
    Try {
      val uri = parseURI(string)
      val storageAccount = parseStorageAccount(parseURI(endpoint))
      val hasContainer = uri.getPath().split("/").filter(!_.isEmpty()).headOption.contains(container)
      def hasEndpoint = parseStorageAccount(uri).contains(storageAccount.get)
      if (hasContainer && !storageAccount.isEmpty && hasEndpoint) {
        ValidBlobPath(uri.getPath.replaceFirst("/" + container, ""))
      } else {
        UnparsableBlobPath(new MalformedURLException(invalidBlobPathMessage(container, endpoint)))
      }
    } recover { case t => UnparsableBlobPath(t) } get
  }
}

class BlobPathBuilder(blobTokenGenerator: BlobTokenGenerator, container: String, endpoint: String) extends PathBuilder {

  val credential: AzureSasCredential = new AzureSasCredential(blobTokenGenerator.getAccessToken)
  val fileSystemConfig: Map[String, Object] = Map((AzureFileSystem.AZURE_STORAGE_SAS_TOKEN_CREDENTIAL, credential),
                                                  (AzureFileSystem.AZURE_STORAGE_FILE_STORES, container))

  def retrieveFilesystem(uri: URI): Try[FileSystem] = {
    Try(FileSystems.getFileSystem(uri)) recover {
      // If no filesystem already exists, this will create a new connection, with the provided configs
      case _: FileSystemNotFoundException => FileSystems.newFileSystem(uri, fileSystemConfig.asJava)
    }
  }

  def build(string: String): Try[BlobPath] = {
    validateBlobPath(string, container, endpoint) match {
      case ValidBlobPath(path) => for {
            fileSystem <- retrieveFilesystem(new URI("azb://?endpoint=" + endpoint))
            nioPath <- Try(fileSystem.getPath(path))
            blobPath = BlobPath(nioPath, endpoint, container)
          } yield blobPath
      case UnparsableBlobPath(errorMessage: Throwable) => Failure(errorMessage)
    }
  }

  override def name: String = "Azure Blob Storage"
}

// Add args for container, storage account name
case class BlobPath private[blob](nioPath: NioPath, endpoint: String, container: String) extends Path {
  override protected def newPath(nioPath: NioPath): Path = BlobPath(nioPath, endpoint, container)

  override def pathAsString: String = List(endpoint, container, nioPath.toString()).mkString("/")

  override def pathWithoutScheme: String = parseURI(endpoint).getHost + "/" + container + "/" + nioPath.toString()
}
