package cromwell.filesystems.blob

import com.azure.core.credential.AzureSasCredential
import cromwell.filesystems.blob.BlobPathBuilder
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.Files

class BlobPathBuilderSpec extends AnyFlatSpec with Matchers{

  it should "parse a URI into a path" in {
    val testString = "https://storageAccount.blob.core.windows.net/container/path/to/file"
    val evalPath = "/path/to/file"
    val pathBuilder = new BlobPathBuilder(new AzureSasCredential("test"), "container", "storageAccount")
    pathBuilder.validateBlobPath(testString) match {
      case BlobPathBuilder.ValidBlobPath(path) => path should equal(evalPath)
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => fail(errorMessage)
    }
  }

  it should "build a blob path from a test string and read a file" in {
    val storageAccount = "coaexternalstorage"
    val store = "inputs"
    val evalPath = "/test/inputFile.txt"
    val sas = "{SAS TOKEN HERE}"
    val endpoint = "https://" + storageAccount + ".blob.core.windows.net"
    val testString = endpoint + "/" + store + evalPath
    val blobPath: BlobPath = new BlobPathBuilder(new AzureSasCredential(sas), store, storageAccount) build testString getOrElse fail()
    blobPath.container should equal(store)
    blobPath.endpoint should equal(endpoint)
    blobPath.pathAsString should equal(testString)
    blobPath.pathWithoutScheme should equal("/" + store + evalPath)
    val is = Files.newInputStream(blobPath.nioPath)
    val fileText = (is.readAllBytes.map(_.toChar)).mkString
    fileText should include ("This is my test file!!!! Did it work?")
  }
}
