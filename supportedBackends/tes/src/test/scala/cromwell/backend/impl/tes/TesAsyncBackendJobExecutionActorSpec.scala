package cromwell.backend.impl.tes

import common.mock.MockSugar
import cromwell.core.logging.JobLogger
import cromwell.core.path
import cromwell.filesystems.blob.{BlobFileSystemManager, BlobPath, WSMBlobSasTokenGenerator}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.Path
import scala.util.{Failure, Success, Try}

class TesAsyncBackendJobExecutionActorSpec extends AnyFlatSpec with Matchers with MockSugar {
  behavior of "TesAsyncBackendJobExecutionActor"

  val fullyQualifiedName = "this.name.is.more.than.qualified"
  val workflowName = "mockWorkflow"
  val someBlobUrl = "https://lz813a3d637adefec2c6e88f.blob.core.windows.net/sc-d8143fd8-aa07-446d-9ba0-af72203f1794/nyxp6c/tes-internal/configuration/supported-vm-sizes"
  val someNotBlobUrl = "https://www.google.com/path/to/exile"
  var index = 0

  val blobInput_0 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someBlobUrl),
    path = someBlobUrl,
    `type` = Option("FILE"),
    content = None
  )
  index = index+1

  val blobInput_1 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someBlobUrl),
    path = someBlobUrl,
    `type` = Option("FILE"),
    content = None
  )
  index = index+1

  val notBlobInput_1 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someNotBlobUrl + index),
    path = someNotBlobUrl + index,
    `type` = Option("FILE"),
    content = None
  )
  index = index+1

  val notBlobInput_2 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someNotBlobUrl + index),
    path = someNotBlobUrl + index,
    `type` = Option("FILE"),
    content = None
  )

  // Mock blob path functionality.
  val testWsmEndpoint = "https://wsm.mock.com/endpoint"
  val testWorkspaceId = "e58ed763-928c-4155-0000-fdbaaadc15f3"
  val testContainerResourceId = "e58ed763-928c-4155-1111-fdbaaadc15f3"

  def generateMockBlobPath: BlobPath = {
    val mockCromwellPath: cromwell.core.path.Path = mock[cromwell.core.path.Path]
    mockCromwellPath.normalize() returns mockCromwellPath
    val mockNioPath: path.NioPath = mock[path.NioPath]
    val mockJavaPath: Path = mock[java.nio.file.Path]
    mockNioPath.toAbsolutePath returns mockJavaPath
    mockNioPath.normalize() returns mockNioPath

    val mockBlobPath = mock[BlobPath]
    mockBlobPath.nioPath returns mockNioPath
    mockBlobPath.toAbsolutePath returns mockCromwellPath
    mockBlobPath.md5 returns "MOCK_MD5"

    val mockTokenGenerator: WSMBlobSasTokenGenerator = mock[WSMBlobSasTokenGenerator]
    val mockFsm: BlobFileSystemManager = mock[BlobFileSystemManager]
    mockTokenGenerator.getWSMSasFetchEndpoint(mockBlobPath) returns Try(s"$testWsmEndpoint/api/workspaces/v1/$testWorkspaceId/resources/controlled/azure/storageContainer/$testContainerResourceId/getSasToken")
    mockFsm.blobTokenGenerator returns mockTokenGenerator
    mockBlobPath.getFilesystemManager returns mockFsm

    mockBlobPath
  }

  def mockPathGetter(pathString: String): Try[cromwell.core.path.Path] = {
    val mockBlobPath = generateMockBlobPath
    val mockCromwellPath: cromwell.core.path.Path = mock[cromwell.core.path.Path]
    val foundBlobPath: Success[BlobPath] = Success(mockBlobPath)
    val foundNonBlobPath: Success[cromwell.core.path.Path] = Success(mockCromwellPath)
    if (pathString.equals(blobInput_0.url.get) || pathString.equals(blobInput_1.url.get)) return foundBlobPath
    foundNonBlobPath
  }

  def mockBlobConverter(pathToConvert: Try[cromwell.core.path.Path]): Try[BlobPath] = {
    //using a stubbed md5 rather than matching on type because type matching of mocked types at runtime causes problems
    if (pathToConvert.get.md5.equals("MOCK_MD5")) pathToConvert.asInstanceOf[Try[BlobPath]] else Failure(new Exception("failed"))
  }

  it should "not return sas endpoint when no blob paths are provided" in {
    val mockLogger: JobLogger = mock[JobLogger]
    val emptyInputs: List[Input] = List()
    val bloblessInputs: List[Input] = List(notBlobInput_1, notBlobInput_2)
    TesAsyncBackendJobExecutionActor.determineWSMSasEndpointFromInputs(emptyInputs, mockPathGetter, mockLogger, mockBlobConverter).isFailure shouldBe true
    TesAsyncBackendJobExecutionActor.determineWSMSasEndpointFromInputs(bloblessInputs, mockPathGetter, mockLogger, mockBlobConverter).isFailure shouldBe true
  }

  it should "return a sas endpoint based on inputs when blob paths are provided" in {
    val mockLogger: JobLogger = mock[JobLogger]
    val expected = s"$testWsmEndpoint/api/workspaces/v1/$testWorkspaceId/resources/controlled/azure/storageContainer/$testContainerResourceId/getSasToken"
    val blobInput: List[Input] = List(blobInput_0)
    val blobInputs: List[Input] = List(blobInput_0, blobInput_1)
    val mixedInputs: List[Input] = List(notBlobInput_1, blobInput_0, blobInput_1)
    TesAsyncBackendJobExecutionActor.determineWSMSasEndpointFromInputs(blobInput, mockPathGetter, mockLogger, mockBlobConverter).get shouldEqual expected
    TesAsyncBackendJobExecutionActor.determineWSMSasEndpointFromInputs(blobInputs, mockPathGetter, mockLogger, mockBlobConverter).get shouldEqual expected
    TesAsyncBackendJobExecutionActor.determineWSMSasEndpointFromInputs(mixedInputs, mockPathGetter, mockLogger,  mockBlobConverter).get shouldEqual expected
  }

  it should "contain expected strings in the bash script" in {

    val mockEnvironmentVariableNameFromWom = "mock_env_var_for_storing_sas_token"
    val expectedEndpoint = s"$testWsmEndpoint/api/workspaces/v1/$testWorkspaceId/resources/controlled/azure/storageContainer/$testContainerResourceId/getSasToken"

    val beginSubstring = "### BEGIN ACQUIRE LOCAL SAS TOKEN ###"
    val endSubstring = "### END ACQUIRE LOCAL SAS TOKEN ###"
    val curlCommandSubstring =
      s"""
        |sas_response_json=$$(curl -s \\
        |                    --retry 3 \\
        |                    --retry-delay 2 \\
        |                    -X POST "$expectedEndpoint" \\
        |                    -H "Content-Type: application/json" \\
        |                    -H "accept: */*" \\
        |                    -H "Authorization: Bearer $${BEARER_TOKEN}")
        |""".stripMargin
    val exportCommandSubstring = s"""export $mockEnvironmentVariableNameFromWom=$$(echo "$${sas_response_json}" | jq -r '.token')"""

    val generatedBashScript = TesAsyncBackendJobExecutionActor.generateLocalizedSasScriptPreamble(mockEnvironmentVariableNameFromWom, expectedEndpoint)

    generatedBashScript should include (beginSubstring)
    generatedBashScript should include (endSubstring)
    generatedBashScript should include (curlCommandSubstring)
    generatedBashScript should include (exportCommandSubstring)
  }
}
