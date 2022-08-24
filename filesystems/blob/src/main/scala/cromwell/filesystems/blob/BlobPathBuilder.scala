package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import com.azure.storage.blob.models.BlobStorageException
import com.azure.storage.blob.nio.AzureFileSystem
import com.google.common.net.UrlEscapers
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.blob.BlobPathBuilder._

import java.io.IOException
import java.net.{MalformedURLException, URI}
import java.nio.file._
import scala.jdk.CollectionConverters._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

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
      def hasEndpoint = storageAccount.map(parseStorageAccount(uri).contains(_)).getOrElse(false)
      if (hasContainer && !storageAccount.isEmpty && hasEndpoint) {
        ValidBlobPath(uri.getPath.replaceFirst("/" + container, ""))
      } else {
        UnparsableBlobPath(new MalformedURLException(invalidBlobPathMessage(container, endpoint)))
      }
    } recover { case t => UnparsableBlobPath(t) } get
  }
}

class BlobPathBuilder(blobTokenGenerator: BlobTokenGenerator, container: String, endpoint: String) extends PathBuilder {

  def build(string: String): Try[BlobPath] = {
    validateBlobPath(string, container, endpoint) match {
      case ValidBlobPath(path) => Try(BlobPath(path, endpoint, container, blobTokenGenerator))
      case UnparsableBlobPath(errorMessage: Throwable) => Failure(errorMessage)
    }
  }
  override def name: String = "Azure Blob Storage"
}

object BlobPath {

  def buildURI(endpoint: String) = new URI("azb://?endpoint=" + endpoint)

  def buildConfigMap(credential: AzureSasCredential, container: String): Map[String, Object] = {
    Map((AzureFileSystem.AZURE_STORAGE_SAS_TOKEN_CREDENTIAL, credential),
      (AzureFileSystem.AZURE_STORAGE_FILE_STORES, container),
      (AzureFileSystem.AZURE_STORAGE_SKIP_INITIAL_CONTAINER_CHECK, java.lang.Boolean.TRUE))
  }

  def findNioPath(path: String, endpoint: String, container: String, blobTokenGenerator: BlobTokenGenerator, attempted: Boolean = false): NioPath = (for {
    fileSystem <- retrieveFilesystem(buildURI(endpoint), container, blobTokenGenerator)
    nioPath <- Try(retrieveFilePath(fileSystem, path))
  } yield nioPath) match {
    case Success(value) => value
    case Failure(exception: IOException) => exception.getCause() match {
      // Azure NIO library wraps blobStorageExceptions in IOExceptions.
      case cause: BlobStorageException => attempted match {
        // If a restart of the filesystem was already attempted, throw the exception that the IO is wrapping
        case true => throw cause
        // This exception indicated that the filesystem was opened successfully, but something is wrong
        // Try closing the filesystem, and opening with a fresh token
        case false => {
          work
          closeFileSystem(buildURI(endpoint))
          findNioPath(path, endpoint, container, blobTokenGenerator, true)
        }
      }
      case _ => throw exception
    }
    case Failure(exception) => throw exception
  }

  def retrieveFilesystem(uri: URI, container: String, blobTokenGenerator: BlobTokenGenerator): Try[FileSystem] = {
    Try(FileSystems.getFileSystem(uri)) recover {
      // If no filesystem already exists, this will create a new connection, with the provided configs
      case _: FileSystemNotFoundException => {
        val fileSystemConfig = buildConfigMap(blobTokenGenerator.getAccessToken, container)
        FileSystems.newFileSystem(uri, fileSystemConfig.asJava)
      }
    }
  }

  def closeFileSystem(uri: URI) = Try(FileSystems.getFileSystem(uri)).map(_.close)

  def retrieveFilePath(fileSystem: FileSystem, path: String): NioPath = {
    val blobPath = fileSystem.getPath(path)
    fileSystem.provider().checkAccess(blobPath)
    blobPath
  }
}
case class BlobPath private[blob](pathString: String, endpoint: String, container: String, blobTokenGenerator: BlobTokenGenerator) extends Path {
  //var token = blobTokenGenerator.getAccessToken
  //var expiry = token.getSignature.split("&").filter(_.startsWith("se")).headOption.map(_.replaceFirst("se=",""))
  override def nioPath: NioPath = BlobPath.findNioPath(path = pathString, endpoint, container, blobTokenGenerator)

  override protected def newPath(nioPath: NioPath): Path = BlobPath(nioPath.toString(), endpoint, container, blobTokenGenerator)

  override def pathAsString: String = List(endpoint, container, nioPath.toString()).mkString("/")

  override def pathWithoutScheme: String = parseURI(endpoint).getHost + "/" + container + "/" + nioPath.toString()
}
