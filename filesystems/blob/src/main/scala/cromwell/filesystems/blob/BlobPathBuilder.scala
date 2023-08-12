package cromwell.filesystems.blob

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

  sealed trait BlobPathValidation
  case class ValidBlobPath(path: String) extends BlobPathValidation
  case class UnparsableBlobPath(errorMessage: Throwable) extends BlobPathValidation

  def invalidBlobPathMessage(container: BlobContainerName, endpoint: EndpointURL) = s"Malformed Blob URL for this builder. Expecting a URL for a container $container and endpoint $endpoint"
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
  def validateBlobPath(string: String, container: BlobContainerName, endpoint: EndpointURL): BlobPathValidation = {
    val blobValidation = for {
      testUri <- parseURI(string)
      endpointUri <- parseURI(endpoint.value)
      testStorageAccount <- parseStorageAccount(testUri)
      endpointStorageAccount <- parseStorageAccount(endpointUri)
      hasContainer = testUri.getPath.split("/").find(_.nonEmpty).contains(container.value)
      hasEndpoint = testStorageAccount.equals(endpointStorageAccount)
      blobPathValidation = (hasContainer && hasEndpoint) match {
        case true => ValidBlobPath(testUri.getPath.replaceFirst("/" + container, ""))
        case false => UnparsableBlobPath(new MalformedURLException(invalidBlobPathMessage(container, endpoint)))
      }
    } yield blobPathValidation
    blobValidation recover { case t => UnparsableBlobPath(t) } get
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
  // Instead, some uploaders such as TES [2] may optionally calculate the MD5 themselves and store it under this key in metadata.
  // N.B. most if not virtually all large files in the wild will NOT have this key populated because they were not created by TES. [3]
  //
  // [0] https://learn.microsoft.com/en-us/azure/storage/blobs/scalability-targets
  // [1] https://learn.microsoft.com/en-us/rest/api/storageservices/version-2019-12-12
  // [2] https://github.com/microsoft/ga4gh-tes/pull/236
  // [3] As of 2023-08 there are zero search engine results for `md5_hashlist_root_hash` and the only sure-thing client is TES
  private val largeBlobFileMetadataKey = "md5_hashlist_root_hash"

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
    fileSystem <- fsm.retrieveFilesystem()
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
        case (None, Some(metadataMd5)) => Option(hexString(metadataMd5.getBytes))
        case (Some(headerMd5Bytes), None) if headerMd5Bytes.isEmpty => None
        // (Some, Some) could happen if an uploader redundantly populates an md5 for a small file.
        // Doesn't seem like an erroneous condition so just choose the native one.
        case (Some(headerMd5Bytes), _) => Option(hexString(headerMd5Bytes))
      }
    }
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
