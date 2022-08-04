package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import cromwell.filesystems.blob.BlobPathBuilder
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.Files
import scala.util.{Failure, Success, Try}

object BlobPathBuilderSpec {
  def buildEndpoint(storageAccount: String) = s"https://$storageAccount.blob.core.windows.net"
}

class BlobPathBuilderSpec extends AnyFlatSpec with Matchers{

  it should "parse a URI into a path" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = "container"
    val evalPath = "/path/to/file"
    val testString = endpoint + "/" + container + evalPath
    BlobPathBuilder.validateBlobPath(testString, container, endpoint) match {
      case BlobPathBuilder.ValidBlobPath(path) => path should equal(evalPath)
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => fail(errorMessage)
    }
  }

  it should "bad storage account fails causes URI to fail parse into a path" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = "container"
    val evalPath = "/path/to/file"
    val testString = BlobPathBuilderSpec.buildEndpoint("badStorageAccount") + container + evalPath
    BlobPathBuilder.validateBlobPath(testString, container, endpoint) match {
      case BlobPathBuilder.ValidBlobPath(path) => fail(s"Valid path: $path found when verifying mismatched storage account")
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => errorMessage.getMessage() should equal(BlobPathBuilder.invalidBlobPathMessage(container, endpoint))
    }
  }

  it should "bad container fails causes URI to fail parse into a path" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = "container"
    val evalPath = "/path/to/file"
    val testString = endpoint + "badContainer" + evalPath
    BlobPathBuilder.validateBlobPath(testString, container, endpoint) match {
      case BlobPathBuilder.ValidBlobPath(path) => fail(s"Valid path: $path found when verifying mismatched container")
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => errorMessage.getMessage() should equal(BlobPathBuilder.invalidBlobPathMessage(container, endpoint))
    }
  }

  it should "build a blob path from a test string and read a file" in {
    val endpoint = "https://saf3d608c45a3ebef88c7c.blob.core.windows.net" //BlobPathBuilderSpec.buildEndpoint("saf3d608c45a3ebef88c7c")
    val endpointHost = "saf3d608c45a3ebef88c7c" //BlobPathBuilder.parseURI(endpoint).getHost
    val store = "sc-f3d608c4-0a0d-4f23-a248-5a3ebef88c7c"
    val evalPath = "/test.txt"
    val sas = "?sv=2021-06-08&srt=sco&spr=https&st=2022-08-03T20%3A48%3A48Z&se=2022-08-03T22%3A03%3A48Z&sr=c&sp=racwdl&sig=r59lmUcW%2BqBaugZA8lvjbr%2Fy2dvvn1hBNgGpYp%2FgQbg%3D"
    val testString = "https://saf3d608c45a3ebef88c7c.blob.core.windows.net/sc-f3d608c4-0a0d-4f23-a248-5a3ebef88c7c/test.txt"
    val blobPathTry: Try[BlobPath] = new BlobPathBuilder(new AzureSasCredential(sas), store, endpoint) build testString
    blobPathTry match {
      case Success(blobPath) => {
        blobPath.container should equal(store)
        blobPath.endpoint should equal(endpoint)
        blobPath.pathAsString should equal(testString)
        blobPath.pathWithoutScheme should equal(endpointHost + "/" + store + evalPath)
        val is = Files.newInputStream(blobPath.nioPath)
        val fileText = (is.readAllBytes.map(_.toChar)).mkString
        fileText should include ("This is a file!")
      }
      case Failure(exception) => {
        fail(exception)
      }
    }
  }
}
