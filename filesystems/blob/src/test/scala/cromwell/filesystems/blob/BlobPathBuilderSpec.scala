package cromwell.filesystems.blob
import common.mock.MockSugar
import org.mockito.Mockito.when
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.util.UUID
import scala.util.{Failure, Try}

object BlobPathBuilderSpec {
  def buildEndpoint(storageAccount: String) = EndpointURL(s"https://$storageAccount.blob.core.windows.net")
}

class BlobPathBuilderSpec extends AnyFlatSpec with Matchers with MockSugar {
  // ValidateBlobPath
  it should "parse a URI into a path" in {
    val endpointInput = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val containerInput = BlobContainerName("container")
    val evalPath = "/path/to/file"
    val testString = endpointInput.value + "/" + containerInput + evalPath
    BlobPathBuilder.validateBlobPath(testString) match {
      case BlobPathBuilder.ValidBlobPath(path, container, endpoint) => {
        path should equal(evalPath)
        container should equal(containerInput)
        endpoint should equal(endpointInput)
      }
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => fail(errorMessage)
    }
  }

  it should "fail to parse an invalid URI into a path" in {
    val endpointInput = EndpointURL("https://storageAccount.bad.host")
    val containerInput = BlobContainerName("container")
    val evalPath = "/path/to/file"
    val testString = endpointInput.value + "/"+ containerInput + evalPath
    BlobPathBuilder.validateBlobPath(testString) match {
      case BlobPathBuilder.ValidBlobPath(_, _, _) => fail("Invalid blob path found as valid")
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => {
        errorMessage.getMessage() should equal(BlobPathBuilder.invalidBlobHostMessage(endpointInput))
      }
    }
  }

  it should "provide a readable error when getting an illegal nioPath" in {
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("container")
    val evalPath = "/path/to/file"
    val exception = new Exception("Failed to do the thing")
    val fsm = mock[BlobFileSystemManager]
    when(fsm.retrieveFilesystem(endpoint, container)).thenReturn(Failure(exception))
    val path = BlobPath(evalPath, endpoint, container)(fsm)
    val testException = Try(path.nioPath).failed.toOption
    testException should contain(exception)
  }

  private def testBlobNioStringCleaning(input: String, expected: String) =
    BlobPath.cleanedNioPathString(input) shouldBe expected

  it should "clean the NIO path string when it has a garbled http protocol" in {
    testBlobNioStringCleaning(
      "https:/lz43.blob.core.windows.net/sc-ebda3e/workspace-services/cbas/terra-app-4628d0e1/test_all_engine_functions/4bb6a0a2-3b07/call-run_read_string/execution/stdout",
      "/workspace-services/cbas/terra-app-4628d0e1/test_all_engine_functions/4bb6a0a2-3b07/call-run_read_string/execution/stdout"
    )
  }

  it should "clean the NIO path string when it has a container name with colon prefix" in {
    testBlobNioStringCleaning(
      "sc-ebda3e:/workspace-services/cbas/terra-app-4628d0e1/test_all_engine_functions/4bb6a0a2-3b07/call-run_read_string/execution/stdout",
      "/workspace-services/cbas/terra-app-4628d0e1/test_all_engine_functions/4bb6a0a2-3b07/call-run_read_string/execution/stdout"
    )
  }

  it should "clean the NIO path string when it's an in-container absolute path" in {
    testBlobNioStringCleaning(
      "/workspace-services/cbas/terra-app-4628d0e1/test_all_engine_functions/4bb6a0a2-3b07/call-run_read_string/execution/stdout",
      "/workspace-services/cbas/terra-app-4628d0e1/test_all_engine_functions/4bb6a0a2-3b07/call-run_read_string/execution/stdout"
    )
  }

  it should "clean the NIO path string when it's the root directory only" in {
    testBlobNioStringCleaning(
      "sc-ebda3e:",
      ""
    )
  }

  //// The below tests are IGNORED because they depend on Azure auth information being present in the environment ////
  private val subscriptionId: SubscriptionId = SubscriptionId(UUID.fromString("62b22893-6bc1-46d9-8a90-806bb3cce3c9"))
  private val endpoint: EndpointURL = BlobPathBuilderSpec.buildEndpoint("coaexternalstorage")
  private val store: BlobContainerName = BlobContainerName("inputs")

  def makeBlobPathBuilder(blobEndpoint: EndpointURL, container: BlobContainerName): BlobPathBuilder = {
    val blobTokenGenerator = NativeBlobSasTokenGenerator(container, blobEndpoint, Some(subscriptionId))
    val fsm = new BlobFileSystemManager(container, blobEndpoint, 10, blobTokenGenerator)
    new BlobPathBuilder(store, endpoint)(fsm)
  }

  ignore should "resolve an absolute path string correctly to a path" in {
    val builder = makeBlobPathBuilder(endpoint, store)
    val rootString = s"${endpoint.value}/${store.value}/cromwell-execution"
    val blobRoot: BlobPath = builder build rootString getOrElse fail()
    blobRoot.toAbsolutePath.pathAsString should equal ("https://coaexternalstorage.blob.core.windows.net/inputs/cromwell-execution")
    val otherFile = blobRoot.resolve("https://coaexternalstorage.blob.core.windows.net/inputs/cromwell-execution/test/inputFile.txt")
    otherFile.toAbsolutePath.pathAsString should equal ("https://coaexternalstorage.blob.core.windows.net/inputs/cromwell-execution/test/inputFile.txt")
  }

  ignore should "build a blob path from a test string and read a file" in {
    val builder = makeBlobPathBuilder(endpoint, store)
    val endpointHost = BlobPathBuilder.parseURI(endpoint.value).map(_.getHost).getOrElse(fail("Could not parse URI"))
    val evalPath = "/test/inputFile.txt"
    val testString = endpoint.value + "/" + store + evalPath
    val blobPath: BlobPath = builder build testString getOrElse fail()

    blobPath.container should equal(store)
    blobPath.endpoint should equal(endpoint)
    blobPath.pathAsString should equal(testString)
    blobPath.pathWithoutScheme should equal(endpointHost + "/" + store + evalPath)
    val is = blobPath.newInputStream()
    val fileText = (is.readAllBytes.map(_.toChar)).mkString
    fileText should include ("This is my test file!!!! Did it work?")
  }

  ignore should "build duplicate blob paths in the same filesystem" in {
    val builder = makeBlobPathBuilder(endpoint, store)
    val evalPath = "/test/inputFile.txt"
    val testString = endpoint.value + "/" + store + evalPath
    val blobPath1: BlobPath = builder build testString getOrElse fail()
    blobPath1.nioPath.getFileSystem.close()
    val blobPath2: BlobPath = builder build testString getOrElse fail()
    blobPath1 should equal(blobPath2)
    val is = blobPath1.newInputStream()
    val fileText = (is.readAllBytes.map(_.toChar)).mkString
    fileText should include ("This is my test file!!!! Did it work?")
  }

  ignore should "resolve a path without duplicating container name" in {
    val builder = makeBlobPathBuilder(endpoint, store)
    val rootString = s"${endpoint.value}/${store.value}/cromwell-execution"
    val blobRoot: BlobPath = builder build rootString getOrElse fail()
    blobRoot.toAbsolutePath.pathAsString should equal ("https://coaexternalstorage.blob.core.windows.net/inputs/cromwell-execution")
    val otherFile = blobRoot.resolve("test/inputFile.txt")
    otherFile.toAbsolutePath.pathAsString should equal ("https://coaexternalstorage.blob.core.windows.net/inputs/cromwell-execution/test/inputFile.txt")
  }

  ignore should "correctly remove a prefix from the blob path" in {
    val builder = makeBlobPathBuilder(endpoint, store)
    val rootString = s"${endpoint.value}/${store.value}/cromwell-execution/"
    val execDirString = s"${endpoint.value}/${store.value}/cromwell-execution/abc123/myworkflow/task1/def4356/execution/"
    val fileString = s"${endpoint.value}/${store.value}/cromwell-execution/abc123/myworkflow/task1/def4356/execution/stdout"
    val blobRoot: BlobPath = builder build rootString getOrElse fail()
    val execDir: BlobPath = builder build execDirString getOrElse fail()
    val blobFile: BlobPath = builder build fileString getOrElse fail()
    blobFile.pathStringWithoutPrefix(blobRoot) should equal ("abc123/myworkflow/task1/def4356/execution/stdout")
    blobFile.pathStringWithoutPrefix(execDir) should equal ("stdout")
    blobFile.pathStringWithoutPrefix(blobFile) should equal ("")
  }

  ignore should "not change a path if it doesn't start with a prefix" in {
    val builder = makeBlobPathBuilder(endpoint, store)
    val otherRootString = s"${endpoint.value}/${store.value}/foobar/"
    val fileString = s"${endpoint.value}/${store.value}/cromwell-execution/abc123/myworkflow/task1/def4356/execution/stdout"
    val otherBlobRoot: BlobPath = builder build otherRootString getOrElse fail()
    val blobFile: BlobPath = builder build fileString getOrElse fail()
    blobFile.pathStringWithoutPrefix(otherBlobRoot) should equal ("/cromwell-execution/abc123/myworkflow/task1/def4356/execution/stdout")
  }
}
