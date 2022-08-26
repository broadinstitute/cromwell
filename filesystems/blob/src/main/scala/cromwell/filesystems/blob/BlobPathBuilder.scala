package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import com.google.common.net.UrlEscapers
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.blob.BlobPathBuilder._

import java.net.{MalformedURLException, URI}
import java.time.Instant
import java.time.temporal.TemporalAmount
import scala.language.postfixOps
import scala.util.{Failure, Try}

object BlobPathBuilder {

  sealed trait BlobPathValidation
  case class ValidBlobPath(path: String) extends BlobPathValidation
  case class UnparsableBlobPath(errorMessage: Throwable) extends BlobPathValidation

  def invalidBlobPathMessage(container: String, endpoint: String) = s"Malformed Blob URL for this builder. Expecting a URL for a container $container and endpoint $endpoint"
  def parseURI(string: String): URI = URI.create(UrlEscapers.urlFragmentEscaper().escape(string))
  def parseStorageAccount(uri: URI): Option[String] = uri.getHost.split("\\.").find(_.nonEmpty)

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
      val hasContainer = uri.getPath.split("/").find(_.nonEmpty).contains(container)
      val hasEndpoint = storageAccount.exists(parseStorageAccount(uri).contains(_))
      if (hasContainer && storageAccount.isDefined && hasEndpoint) {
        ValidBlobPath(uri.getPath.replaceFirst("/" + container, ""))
      } else {
        UnparsableBlobPath(new MalformedURLException(invalidBlobPathMessage(container, endpoint)))
      }
    } recover { case t => UnparsableBlobPath(t) } get
  }
}

class BlobPathBuilder(fsm: FileSystemManager, container: String, endpoint: String) extends PathBuilder {

  def build(string: String): Try[BlobPath] = {
    validateBlobPath(string, container, endpoint) match {
      case ValidBlobPath(path) => Try(BlobPath(path, endpoint, container, fsm))
      case UnparsableBlobPath(errorMessage: Throwable) => Failure(errorMessage)
    }
  }
  override def name: String = "Azure Blob Storage"
}

case class BlobPath private[blob](pathString: String, endpoint: String, container: String, fsm: FileSystemManager) extends Path {
  //var token = blobTokenGenerator.getAccessToken
  //var expiry = token.getSignature.split("&").filter(_.startsWith("se")).headOption.map(_.replaceFirst("se=",""))
  override def nioPath: NioPath = findNioPath(path = pathString, endpoint, container)

  override protected def newPath(nioPath: NioPath): Path = BlobPath(nioPath.toString, endpoint, container, fsm)

  override def pathAsString: String = List(endpoint, container, nioPath.toString).mkString("/")

  override def pathWithoutScheme: String = parseURI(endpoint).getHost + "/" + container + "/" + nioPath.toString

  def findNioPath(path: String, endpoint: String, container: String): NioPath = (for {
    fileSystem <- fsm.retrieveFilesystem()
    nioPath = fileSystem.getPath(path)
  } yield nioPath).get
}

case class TokenExpiration(token: AzureSasCredential, buffer: TemporalAmount) {
  val expiry = for {
    expiryString <- token.getSignature.split("&").find(_.startsWith("se")).map(_.replaceFirst("se=","")).map(_.replace("%3A", ":"))
    instant = Instant.parse(expiryString)
  } yield instant

  def hasTokenExpired: Boolean = expiry.exists(_.isAfter(Instant.now.plus(buffer)))
}
