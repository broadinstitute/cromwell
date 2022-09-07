package cromwell.filesystems.blob

import com.google.common.net.UrlEscapers
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.blob.BlobPathBuilder._

import java.net.{MalformedURLException, URI}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object BlobPathBuilder {

  sealed trait BlobPathValidation
  case class ValidBlobPath(path: String) extends BlobPathValidation
  case class UnparsableBlobPath(errorMessage: Throwable) extends BlobPathValidation

  def invalidBlobPathMessage(container: BlobContainerName, endpoint: EndpointURL) = s"Malformed Blob URL for this builder. Expecting a URL for a container $container and endpoint $endpoint"
  def parseURI(string: String): URI = URI.create(UrlEscapers.urlFragmentEscaper().escape(string))
  def parseStorageAccount(uri: URI): Try[StorageAccountName] = uri.getHost.split("\\.").find(_.nonEmpty).map(StorageAccountName(_)).fold[Try[StorageAccountName]](Failure(new Exception("Could not parse storage account")))(Success(_))

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
  def validateBlobPath(string: String, container: BlobContainerName, endpoint: EndpointURL): BlobPathValidation = {
    Try {
      val uri = parseURI(string)
      val storageAccount = parseStorageAccount(parseURI(endpoint.value))
      val hasContainer = uri.getPath.split("/").find(_.nonEmpty).contains(container.value)
      val hasEndpoint = storageAccount.toOption.exists(parseStorageAccount(uri).toOption.contains(_))
      if (hasContainer && storageAccount.isSuccess && hasEndpoint) {
        ValidBlobPath(uri.getPath.replaceFirst("/" + container, ""))
      } else {
        UnparsableBlobPath(new MalformedURLException(invalidBlobPathMessage(container, endpoint)))
      }
    } recover { case t => UnparsableBlobPath(t) } get
  }
}

class BlobPathBuilder(container: BlobContainerName, endpoint: EndpointURL)(private val fsm: BlobFileSystemManager) extends PathBuilder {

  def build(string: String): Try[BlobPath] = {
    validateBlobPath(string, container, endpoint) match {
      case ValidBlobPath(path) => Try(BlobPath(path, endpoint, container)(fsm))
      case UnparsableBlobPath(errorMessage: Throwable) => Failure(errorMessage)
    }
  }
  override def name: String = "Azure Blob Storage"
}

case class BlobPath private[blob](pathString: String, endpoint: EndpointURL, container: BlobContainerName)(private val fsm: BlobFileSystemManager) extends Path {
  override def nioPath: NioPath = findNioPath(pathString)

  override protected def newPath(nioPath: NioPath): Path = BlobPath(nioPath.toString, endpoint, container)(fsm)

  override def pathAsString: String = List(endpoint, container, nioPath.toString).mkString("/")

  override def pathWithoutScheme: String = parseURI(endpoint.value).getHost + "/" + container + "/" + nioPath.toString

  private def findNioPath(path: String): NioPath = (for {
    fileSystem <- fsm.retrieveFilesystem()
    nioPath = fileSystem.getPath(path)
  // This is purposefully an unprotected get because the NIO API needing an unwrapped path object.
  // If an error occurs the api expects a thrown exception
  } yield nioPath).get
}
