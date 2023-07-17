package cromwell.filesystems.blob

import com.azure.storage.blob.nio.AzureBlobFileAttributes
import com.google.common.net.UrlEscapers
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.blob.BlobPathBuilder._

import java.net.{MalformedURLException, URI}
import java.nio.file.{Files, LinkOption}
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object BlobPathBuilder {
  val HOSTNAME_SUFFIX = ".blob.core.windows.net"
  sealed trait BlobPathValidation
  case class ValidBlobPath(path: String, container: BlobContainerName, endpoint: EndpointURL) extends BlobPathValidation
  case class UnparsableBlobPath(errorMessage: Throwable) extends BlobPathValidation

  def invalidBlobHostMessage(endpoint: EndpointURL) = s"Malformed Blob URL for this builder: The endpoint $endpoint doesn't contain the expected host string '{SA}.blob.core.windows.net/'"
  def invalidBlobContainerMessage(endpoint: EndpointURL) = s"Malformed Blob URL for this builder: Could not parse container"
  def parseURI(string: String): Try[URI] = Try(URI.create(UrlEscapers.urlFragmentEscaper().escape(string)))
  def parseStorageAccount(uri: URI): Try[StorageAccountName] = uri.getHost.split("\\.").find(_.nonEmpty).map(StorageAccountName(_))
      .map(Success(_)).getOrElse(Failure(new Exception("Could not parse storage account")))

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
 def validateBlobPath(string: String): BlobPathValidation = {
    val blobValidation = for {
      testUri <- parseURI(string)
      testEndpoint = EndpointURL(testUri.getScheme + "://" + testUri.getHost())
      testStorageAccount <- parseStorageAccount(testUri)
      testContainer = testUri.getPath.split("/").find(_.nonEmpty)
      isBlobHost = testUri.getHost().contains(HOSTNAME_SUFFIX) && testUri.getScheme().contains("https")
      blobPathValidation = (isBlobHost, testContainer) match {
        case (true, Some(container)) => ValidBlobPath(
            testUri.getPath.replaceFirst("/" + container, ""),
            BlobContainerName(container),
            testEndpoint)
        case (false, _) => UnparsableBlobPath(new MalformedURLException(invalidBlobHostMessage(testEndpoint)))
        case (true, None) => UnparsableBlobPath(new MalformedURLException(invalidBlobContainerMessage(testEndpoint)))
      }
    } yield blobPathValidation
    blobValidation recover { case t => UnparsableBlobPath(t) } get
  }
}

class BlobPathBuilder()(private val fsm: BlobFileSystemManager) extends PathBuilder {

  def build(string: String): Try[BlobPath] = {
    validateBlobPath(string) match {
      case ValidBlobPath(path, container, endpoint) => Try(BlobPath(path, endpoint, container)(fsm))
      case UnparsableBlobPath(errorMessage: Throwable) => Failure(errorMessage)
    }
  }
  override def name: String = "Azure Blob Storage"
}

object BlobPath {
  // The Azure NIO library uses `{containerName}:` as the root of the path (treating the blob container within
  // the storage account similarly to a drive within a computer). This doesn't work well for our need to easily
  // transfer back and forth to and from the blob URL format. It also causes the library to garble full http://
  // paths that it receives (it interprets `http` as the container name); it transforms them to http:/<remainder of path>
  //
  // We transform these library-generated paths in two steps:
  // 1) If the path starts with http:/ (single slash!) transform it to the containerName:<path inside container>
  //    format the library expects
  // 2) If the path looks like <container>:<path>, strip off the <container>: to leave the absolute path inside the container.
  private val brokenPathRegex = "https:/([a-z0-9]+).blob.core.windows.net/([-a-zA-Z0-9]+)/(.*)".r
  def cleanedNioPathString(nioString: String): String = {
    val pathStr = nioString match {
      case brokenPathRegex(_, containerName, pathInContainer) =>
        s"${containerName}:/${pathInContainer}"
      case _ => nioString
    }
    pathStr.substring(pathStr.indexOf(":")+1)
  }

  def apply(nioPath: NioPath,
            endpoint: EndpointURL,
            container: BlobContainerName,
            fsm: BlobFileSystemManager): BlobPath = {
    BlobPath(cleanedNioPathString(nioPath.toString), endpoint, container)(fsm)
  }
}

case class BlobPath private[blob](pathString: String, endpoint: EndpointURL, container: BlobContainerName)(private val fsm: BlobFileSystemManager) extends Path {
  override def nioPath: NioPath = findNioPath(pathString)

  override protected def newPath(nioPath: NioPath): Path = BlobPath(nioPath, endpoint, container, fsm)

  override def pathAsString: String = List(endpoint, container, pathString.stripPrefix("/")).mkString("/")

  //This is purposefully an unprotected get because if the endpoint cannot be parsed this should fail loudly rather than quietly
  override def pathWithoutScheme: String = parseURI(endpoint.value).map(u => List(u.getHost, container, pathString.stripPrefix("/")).mkString("/")).get

  private def findNioPath(path: String): NioPath = (for {
    fileSystem <- fsm.retrieveFilesystem(endpoint, container)
    // The Azure NIO library uses `{container}:` to represent the root of the path
    nioPath = fileSystem.getPath(s"${container.value}:", path)
  // This is purposefully an unprotected get because the NIO API needing an unwrapped path object.
  // If an error occurs the api expects a thrown exception
  } yield nioPath).get

  def blobFileAttributes: Try[AzureBlobFileAttributes] =
    Try(Files.readAttributes(nioPath, classOf[AzureBlobFileAttributes]))

  def md5HexString: Try[Option[String]] = {
    blobFileAttributes.map(h =>
      Option(h.blobHttpHeaders().getContentMd5) match {
        case None => None
        case Some(arr) if arr.isEmpty => None
        // Convert the bytes to a hex-encoded string. Note that this value
        // is rendered in base64 in the Azure web portal.
        case Some(bytes) => Option(bytes.map("%02x".format(_)).mkString)
      }
    )
  }

  /**
    * Return the pathString of this BlobPath, with the given prefix removed if this path shares that
    * prefix.
    */
  def pathStringWithoutPrefix(prefix: Path): String = {
    if (this.startsWith(prefix)) {
      prefix.relativize(this) match {
        case b: BlobPath => b.pathString // path inside the blob container
        case p: Path => p.pathAsString // full path
      }
    }
    else pathString
  }

  override def getSymlinkSafePath(options: LinkOption*): Path  = toAbsolutePath
}
