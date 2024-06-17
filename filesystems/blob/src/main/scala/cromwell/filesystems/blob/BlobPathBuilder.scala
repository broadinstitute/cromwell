package cromwell.filesystems.blob

import akka.http.scaladsl.model.Uri
import com.azure.storage.blob.nio.AzureBlobFileAttributes
import com.google.common.net.UrlEscapers
import cromwell.core.path.{NioPath, Path, PathBuilder}
import cromwell.filesystems.blob.BlobPathBuilder._

import java.net.{MalformedURLException, URI}
import java.nio.file.{Files, LinkOption}
import scala.jdk.CollectionConverters._
import scala.language.postfixOps
import scala.util.{Failure, Success, Try}

object BlobPathBuilder {
  private val blobHostnameSuffix = ".blob.core.windows.net"
  sealed trait BlobPathValidation
  case class ValidBlobPath(path: String, container: BlobContainerName, endpoint: EndpointURL) extends BlobPathValidation
  case class UnparsableBlobPath(errorMessage: Throwable) extends BlobPathValidation

  def invalidBlobHostMessage(endpoint: EndpointURL) =
    s"Malformed Blob URL for this builder: The endpoint $endpoint doesn't contain the expected host string '{SA}.blob.core.windows.net/'"
  def invalidBlobContainerMessage(endpoint: EndpointURL) =
    s"Malformed Blob URL for this builder: Could not parse container"
  val externalToken =
    "Rejecting pre-signed SAS URL so that filesystem selection falls through to HTTP filesystem"
  def parseURI(string: String): Try[URI] = Try(URI.create(UrlEscapers.urlFragmentEscaper().escape(string)))
  def parseStorageAccount(uri: URI): Try[StorageAccountName] = uri.getHost
    .split("\\.")
    .find(_.nonEmpty)
    .map(StorageAccountName(_))
    .map(Success(_))
    .getOrElse(Failure(new Exception("Could not parse storage account")))

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
      testEndpoint = EndpointURL(testUri.getScheme + "://" + testUri.getHost)
      _ <- parseStorageAccount(testUri)
      testContainer = testUri.getPath.split("/").find(_.nonEmpty)
      isBlobHost = testUri.getHost.contains(blobHostnameSuffix) && testUri.getScheme.contains("https")
      hasToken = hasSasToken(string)
      blobPathValidation = (isBlobHost, testContainer, hasToken) match {
        case (true, Some(container), false) =>
          ValidBlobPath(testUri.getPath.replaceFirst("/" + container, ""), BlobContainerName(container), testEndpoint)
        case (false, _, false) =>
          UnparsableBlobPath(new MalformedURLException(invalidBlobHostMessage(testEndpoint)))
        case (true, None, false) =>
          UnparsableBlobPath(new MalformedURLException(invalidBlobContainerMessage(testEndpoint)))
        case (_, _, true) =>
          UnparsableBlobPath(new IllegalArgumentException(externalToken))
      }
    } yield blobPathValidation
    blobValidation recover { case t => UnparsableBlobPath(t) } get
  }

  private def hasSasToken(uri: Uri) = {
    // These keys are required for all SAS tokens.
    // https://learn.microsoft.com/en-us/rest/api/storageservices/create-service-sas#construct-a-service-sas
    val SignedVersionKey = "sv"
    val SignatureKey = "sig"

    val query = uri.query().toMap
    query.isDefinedAt(SignedVersionKey) && query.isDefinedAt(SignatureKey)
  }
}

class BlobPathBuilder()(private val fsm: BlobFileSystemManager) extends PathBuilder {

  def build(string: String): Try[BlobPath] =
    validateBlobPath(string) match {
      case ValidBlobPath(path, container, endpoint) => Try(BlobPath(path, endpoint, container)(fsm))
      case UnparsableBlobPath(errorMessage: Throwable) => Failure(errorMessage)
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

  // Blob files larger than 5 GB upload in parallel parts [0][1] and do not get a native `CONTENT-MD5` property.
  // Instead, some uploaders such as TES [2] calculate the md5 themselves and store it under this key in metadata.
  // They do this for all files they touch, regardless of size, and the root/metadata property is authoritative over native.
  //
  // N.B. most if not virtually all large files in the wild will NOT have this key populated because they were not created
  // by TES or its associated upload utility [3].
  //
  // [0] https://learn.microsoft.com/en-us/azure/storage/blobs/scalability-targets
  // [1] https://learn.microsoft.com/en-us/rest/api/storageservices/version-2019-12-12
  // [2] https://github.com/microsoft/ga4gh-tes/blob/03feb746bb961b72fa91266a56db845e3b31be27/src/Tes.Runner/Transfer/BlobBlockApiHttpUtils.cs#L25
  // [3] https://github.com/microsoft/ga4gh-tes/blob/main/src/Tes.RunnerCLI/scripts/roothash.sh
  private val largeBlobFileMetadataKey = "md5_4mib_hashlist_root_hash"

  def cleanedNioPathString(nioString: String): String = {
    val pathStr = nioString match {
      case brokenPathRegex(_, containerName, pathInContainer) =>
        s"${containerName}:/${pathInContainer}"
      case _ => nioString
    }
    pathStr.substring(pathStr.indexOf(":") + 1)
  }

  def apply(nioPath: NioPath,
            endpoint: EndpointURL,
            container: BlobContainerName,
            fsm: BlobFileSystemManager
  ): BlobPath =
    BlobPath(cleanedNioPathString(nioPath.toString), endpoint, container)(fsm)
}

case class BlobPath private[blob] (pathString: String, endpoint: EndpointURL, container: BlobContainerName)(
  private val fsm: BlobFileSystemManager
) extends Path {
  override def nioPath: NioPath = findNioPath(pathString)

  override protected def newPath(nioPath: NioPath): Path = BlobPath(nioPath, endpoint, container, fsm)

  override def pathAsString: String = List(endpoint, container, pathString.stripPrefix("/")).mkString("/")

  // This is purposefully an unprotected get because if the endpoint cannot be parsed this should fail loudly rather than quietly
  override def pathWithoutScheme: String =
    parseURI(endpoint.value).map(u => List(u.getHost, container, pathString.stripPrefix("/")).mkString("/")).get

  private def findNioPath(path: String): NioPath = (for {
    fileSystem <- fsm.retrieveFilesystem(endpoint, container)
    // The Azure NIO library uses `{container}:` to represent the root of the path
    nioPath = fileSystem.getPath(s"${container.value}:", path)
    // This is purposefully an unprotected get because the NIO API needing an unwrapped path object.
    // If an error occurs the api expects a thrown exception
  } yield nioPath).get

  def blobFileAttributes: Try[AzureBlobFileAttributes] =
    Try(Files.readAttributes(nioPath, classOf[AzureBlobFileAttributes]))

  def blobFileMetadata: Try[Option[Map[String, String]]] = blobFileAttributes.map { attrs =>
    // `metadata()` has a documented `null` case
    Option(attrs.metadata()).map(_.asScala.toMap)
  }

  def md5HexString: Try[Option[String]] = {
    def md5FromMetadata: Option[String] = (blobFileMetadata map { maybeMetadataMap: Option[Map[String, String]] =>
      maybeMetadataMap flatMap { metadataMap: Map[String, String] =>
        metadataMap.get(BlobPath.largeBlobFileMetadataKey)
      }
    }).toOption.flatten

    // Convert the bytes to a hex-encoded string. Note that the value
    // is rendered in base64 in the Azure web portal.
    def hexString(bytes: Array[Byte]): String = bytes.map("%02x".format(_)).mkString

    blobFileAttributes.map { attr: AzureBlobFileAttributes =>
      (Option(attr.blobHttpHeaders().getContentMd5), md5FromMetadata) match {
        case (None, None) => None
        // (Some, Some) will happen for all <5 GB files uploaded by TES. Per Microsoft 2023-08-15 the
        // root/metadata algorithm emits different values than the native algorithm and we should
        // always choose metadata for consistency with larger files that only have that one.
        case (_, Some(metadataMd5)) => Option(metadataMd5)
        case (Some(headerMd5Bytes), None) if headerMd5Bytes.isEmpty => None
        case (Some(headerMd5Bytes), None) => Option(hexString(headerMd5Bytes))
      }
    }
  }

  /**
    * Return the pathString of this BlobPath, with the given prefix removed if this path shares that
    * prefix.
    */
  def pathStringWithoutPrefix(prefix: Path): String =
    if (this.startsWith(prefix)) {
      prefix.relativize(this) match {
        case b: BlobPath => b.pathString // path inside the blob container
        case p: Path => p.pathAsString // full path
      }
    } else pathString

  /**
    * Returns the path relative to the container root.
    * For example, https://{storageAccountName}.blob.core.windows.net/{containerid}/path/to/my/file
    * will be returned as path/to/my/file.
    * @return Path string relative to the container root.
    */
  def pathWithoutContainer: String = pathString

  def getFilesystemManager: BlobFileSystemManager = fsm

  override def getSymlinkSafePath(options: LinkOption*): Path = toAbsolutePath

}
