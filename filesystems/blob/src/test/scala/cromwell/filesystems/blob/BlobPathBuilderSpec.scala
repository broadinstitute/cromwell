package cromwell.filesystems.blob
import common.mock.MockSugar
import org.mockito.Mockito.when
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import scala.util.{Failure, Try}

object BlobPathBuilderSpec {
  def buildEndpoint(storageAccount: String) = EndpointURL(s"https://$storageAccount.blob.core.windows.net")
}

class BlobPathBuilderSpec extends AnyFlatSpec with Matchers with MockSugar {
  // ValidateBlobPath
  it should "parse a URI into a path" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("container")
    val evalPath = "/path/to/file"
    val testString = endpoint.value + "/" + container + evalPath
    BlobPathBuilder.validateBlobPath(testString, container, endpoint) match {
      case BlobPathBuilder.ValidBlobPath(path) => path should equal(evalPath)
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => fail(errorMessage)
    }
  }

  it should "bad storage account fails causes URI to fail parse into a path" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("container")
    val evalPath = "/path/to/file"
    val testString = BlobPathBuilderSpec.buildEndpoint("badStorageAccount").value + container.value + evalPath
    BlobPathBuilder.validateBlobPath(testString, container, endpoint) match {
      case BlobPathBuilder.ValidBlobPath(path) => fail(s"Valid path: $path found when verifying mismatched storage account")
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => errorMessage.getMessage should equal(BlobPathBuilder.invalidBlobPathMessage(container, endpoint))
    }
  }

  it should "bad container fails causes URI to fail parse into a path" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("container")
    val evalPath = "/path/to/file"
    val testString = endpoint.value + "badContainer" + evalPath
    BlobPathBuilder.validateBlobPath(testString, container, endpoint) match {
      case BlobPathBuilder.ValidBlobPath(path) => fail(s"Valid path: $path found when verifying mismatched container")
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => errorMessage.getMessage should equal(BlobPathBuilder.invalidBlobPathMessage(container, endpoint))
    }
  }

  it should "provide a readable error when getting an illegal nioPath" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("container")
    val evalPath = "/path/to/file"
    val exception = new Exception("Failed to do the thing")
    val fsm = mock[BlobFileSystemManager]
    when(fsm.retrieveFilesystem()).thenReturn(Failure(exception))
    val path = BlobPath(evalPath, endpoint, container)(fsm)
    val testException = Try(path.nioPath).failed.toOption
    testException should contain(exception)
  }

  ignore should "build a blob path from a test string and read a file" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("coaexternalstorage")
    val endpointHost = BlobPathBuilder.parseURI(endpoint.value).map(_.getHost).getOrElse(fail("Could not parse URI"))
    val store = BlobContainerName("inputs")
    val evalPath = "/test/inputFile.txt"
    val blobTokenGenerator = NativeBlobSasTokenGenerator(store, endpoint)
    val fsm: BlobFileSystemManager = new BlobFileSystemManager(store, endpoint, 10L, blobTokenGenerator)
    val testString = endpoint.value + "/" + store + evalPath
    val blobPath: BlobPath = new BlobPathBuilder(store, endpoint)(fsm) build testString getOrElse fail()

    blobPath.container should equal(store)
    blobPath.endpoint should equal(endpoint)
    blobPath.pathAsString should equal(testString)
    blobPath.pathWithoutScheme should equal(endpointHost + "/" + store + evalPath)
    val is = blobPath.newInputStream()
    val fileText = (is.readAllBytes.map(_.toChar)).mkString
    fileText should include ("This is my test file!!!! Did it work?")
  }

  ignore should "build duplicate blob paths in the same filesystem" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("coaexternalstorage")
    val store = BlobContainerName("inputs")
    val evalPath = "/test/inputFile.txt"
    val blobTokenGenerator = NativeBlobSasTokenGenerator(store, endpoint)
    val fsm: BlobFileSystemManager = new BlobFileSystemManager(store, endpoint, 10, blobTokenGenerator)
    val testString = endpoint.value + "/" + store + evalPath
    val blobPath1: BlobPath = new BlobPathBuilder(store, endpoint)(fsm) build testString getOrElse fail()
    blobPath1.nioPath.getFileSystem.close()
    val blobPath2: BlobPath = new BlobPathBuilder(store, endpoint)(fsm) build testString getOrElse fail()
    blobPath1 should equal(blobPath2)
    val is = blobPath1.newInputStream()
    val fileText = (is.readAllBytes.map(_.toChar)).mkString
    fileText should include ("This is my test file!!!! Did it work?")
  }

  ignore should "resolve a path without duplicating container name" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("coaexternalstorage")
    val store = BlobContainerName("inputs")
    val blobTokenGenerator = NativeBlobSasTokenGenerator(store, endpoint)
    val fsm: BlobFileSystemManager = new BlobFileSystemManager(store, endpoint, 10, blobTokenGenerator)

    val rootString = s"${endpoint.value}/${store.value}/cromwell-execution"
    val blobRoot: BlobPath = new BlobPathBuilder(store, endpoint)(fsm) build rootString getOrElse fail()
    blobRoot.toAbsolutePath.pathAsString should equal ("https://coaexternalstorage.blob.core.windows.net/inputs/cromwell-execution")
    val otherFile = blobRoot.resolve("test/inputFile.txt")
    otherFile.toAbsolutePath.pathAsString should equal ("https://coaexternalstorage.blob.core.windows.net/inputs/cromwell-execution/test/inputFile.txt")
  }
}
