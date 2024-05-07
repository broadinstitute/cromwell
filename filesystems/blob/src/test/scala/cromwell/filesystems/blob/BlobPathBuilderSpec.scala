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
    val endpoint = BlobPathBuilderSpec.buildEndpoint("storageAccount")
    val container = BlobContainerName("container")
    val evalPath = "/path/to/file"
    val testString = endpoint.value + "/" + container + evalPath
    BlobPathBuilder.validateBlobPath(testString) match {
      case BlobPathBuilder.ValidBlobPath(path, parsedContainer, parsedEndpoint) =>
        path should equal(evalPath)
        parsedContainer should equal(container)
        parsedEndpoint should equal(endpoint)
      case BlobPathBuilder.UnparsableBlobPath(errorMessage) => fail(errorMessage)
    }
  }

  it should "reject a path that is otherwise valid, but has a preexisting SAS token" in {
    import cromwell.filesystems.blob.BlobPathBuilder.UnparsableBlobPath

    // The `.asInstanceOf[UnparsableBlobPath].errorMessage.getMessage` malarkey is necessary
    // because Java exceptions compare by reference, while strings are by value

    val sasBlob = "https://lz304a1e79fd7359e5327eda.blob.core.windows.net/sc-705b830a-d699-478e-9da6-49661b326e77" +
      "?sv=2021-12-02&spr=https&st=2023-12-13T20%3A27%3A55Z&se=2023-12-14T04%3A42%3A55Z&sr=c&sp=racwdlt&sig=blah&rscd=foo"
    BlobPathBuilder.validateBlobPath(sasBlob).asInstanceOf[UnparsableBlobPath].errorMessage.getMessage should equal(
      UnparsableBlobPath(
        new IllegalArgumentException(
          "Rejecting pre-signed SAS URL so that filesystem selection falls through to HTTP filesystem"
        )
      ).errorMessage.getMessage
    )
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

  // The following tests use the `centaurtesting` account injected into CI. They depend on access to the
  // container specified below. You may need to log in to az cli locally to get them to pass.
  private val subscriptionId: SubscriptionId = SubscriptionId(UUID.fromString("62b22893-6bc1-46d9-8a90-806bb3cce3c9"))

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
  private val endpoint: EndpointURL = BlobPathBuilderSpec.buildEndpoint("centaurtesting")
  private val container: BlobContainerName = BlobContainerName("test-blob")

  def makeBlobPathBuilder(): BlobPathBuilder = {
    val blobTokenGenerator = NativeBlobSasTokenGenerator(Some(subscriptionId))
    val fsm = new BlobFileSystemManager(10, blobTokenGenerator)
    new BlobPathBuilder()(fsm)
  }

  private def testBlobNioStringCleaning(input: String, expected: String) =
    BlobPath.cleanedNioPathString(input) shouldBe expected

  it should "read md5 from small files <5g" in {
    val builder = makeBlobPathBuilder()
    val evalPath = "/testRead.txt"
    val testString = endpoint.value + "/" + container + evalPath
    val blobPath1: BlobPath = (builder build testString).get
    blobPath1.md5HexString.get should equal(Option("31ae06882d06a20e01ba1ac961ce576c"))
  }

  it should "read md5 from large files >5g" in {
    val builder = makeBlobPathBuilder()
    val evalPath = "/Rocky-9.2-aarch64-dvd.iso"
    val testString = endpoint.value + "/" + container + evalPath
    val blobPath1: BlobPath = (builder build testString).get
    blobPath1.md5HexString.toOption.get should equal(Some("13cb09331d2d12c0f476f81c672a4319"))
  }

  it should "choose the root/metadata md5 over the native md5 for files that have both" in {
    val builder = makeBlobPathBuilder()
    val evalPath = "/redundant_md5_test.txt"
    val testString = endpoint.value + "/" + container + evalPath
    val blobPath1: BlobPath = (builder build testString).get
    blobPath1.md5HexString.toOption.get should equal(Some("021c7cc715ec82292bb9b925f9ca44d3"))
  }

  it should "gracefully return `None` when neither hash is found" in {
    val builder = makeBlobPathBuilder()
    val evalPath = "/no_md5_test.txt"
    val testString = endpoint.value + "/" + container + evalPath
    val blobPath1: BlobPath = (builder build testString).get
    blobPath1.md5HexString.get should equal(None)
  }

  it should "resolve an absolute path string correctly to a path" in {
    val builder = makeBlobPathBuilder()
    val rootString = s"${endpoint.value}/${container.value}/cromwell-execution"
    val blobRoot: BlobPath = builder build rootString getOrElse fail()
    blobRoot.toAbsolutePath.pathAsString should equal(
      "https://centaurtesting.blob.core.windows.net/test-blob/cromwell-execution"
    )
    val otherFile =
      blobRoot.resolve("https://centaurtesting.blob.core.windows.net/test-blob/cromwell-execution/test/inputFile.txt")
    otherFile.toAbsolutePath.pathAsString should equal(
      "https://centaurtesting.blob.core.windows.net/test-blob/cromwell-execution/test/inputFile.txt"
    )
  }

  it should "build a blob path from a test string and read a file" in {
    val builder = makeBlobPathBuilder()
    val endpointHost = BlobPathBuilder.parseURI(endpoint.value).map(_.getHost).getOrElse(fail("Could not parse URI"))
    val evalPath = "/test/inputFile.txt"
    val testString = endpoint.value + "/" + container + evalPath
    val blobPath: BlobPath = builder build testString getOrElse fail()

    blobPath.container should equal(container)
    blobPath.endpoint should equal(endpoint)
    blobPath.pathAsString should equal(testString)
    blobPath.pathWithoutScheme should equal(endpointHost + "/" + container + evalPath)
    val is = blobPath.newInputStream()
    val fileText = is.readAllBytes.map(_.toChar).mkString
    fileText should include("This is my test file!!!! Did it work?")
  }

  it should "build duplicate blob paths in the same filesystem" in {
    val builder = makeBlobPathBuilder()
    val evalPath = "/test/inputFile.txt"
    val testString = endpoint.value + "/" + container + evalPath
    val blobPath1: BlobPath = builder build testString getOrElse fail()
    blobPath1.nioPath.getFileSystem.close()
    val blobPath2: BlobPath = builder build testString getOrElse fail()
    blobPath1 should equal(blobPath2)
    val is = blobPath1.newInputStream()
    val fileText = is.readAllBytes.map(_.toChar).mkString
    fileText should include("This is my test file!!!! Did it work?")
  }

  it should "resolve a path without duplicating container name" in {
    val builder = makeBlobPathBuilder()
    val rootString = s"${endpoint.value}/${container.value}/cromwell-execution"
    val blobRoot: BlobPath = builder build rootString getOrElse fail()
    blobRoot.toAbsolutePath.pathAsString should equal(
      "https://centaurtesting.blob.core.windows.net/test-blob/cromwell-execution"
    )
    val otherFile = blobRoot.resolve("test/inputFile.txt")
    otherFile.toAbsolutePath.pathAsString should equal(
      "https://centaurtesting.blob.core.windows.net/test-blob/cromwell-execution/test/inputFile.txt"
    )
  }

  it should "correctly remove a prefix from the blob path" in {
    val builder = makeBlobPathBuilder()
    val rootString = s"${endpoint.value}/${container.value}/cromwell-execution/"
    val execDirString =
      s"${endpoint.value}/${container.value}/cromwell-execution/abc123/myworkflow/task1/def4356/execution/"
    val fileString =
      s"${endpoint.value}/${container.value}/cromwell-execution/abc123/myworkflow/task1/def4356/execution/stdout"
    val blobRoot: BlobPath = builder build rootString getOrElse fail()
    val execDir: BlobPath = builder build execDirString getOrElse fail()
    val blobFile: BlobPath = builder build fileString getOrElse fail()
    blobFile.pathStringWithoutPrefix(blobRoot) should equal("abc123/myworkflow/task1/def4356/execution/stdout")
    blobFile.pathStringWithoutPrefix(execDir) should equal("stdout")
    blobFile.pathStringWithoutPrefix(blobFile) should equal("")
  }

  it should "not change a path if it doesn't start with a prefix" in {
    val builder = makeBlobPathBuilder()
    val otherRootString = s"${endpoint.value}/${container.value}/foobar/"
    val fileString =
      s"${endpoint.value}/${container.value}/cromwell-execution/abc123/myworkflow/task1/def4356/execution/stdout"
    val otherBlobRoot: BlobPath = builder build otherRootString getOrElse fail()
    val blobFile: BlobPath = builder build fileString getOrElse fail()
    blobFile.pathStringWithoutPrefix(otherBlobRoot) should equal(
      "/cromwell-execution/abc123/myworkflow/task1/def4356/execution/stdout"
    )
  }
}
