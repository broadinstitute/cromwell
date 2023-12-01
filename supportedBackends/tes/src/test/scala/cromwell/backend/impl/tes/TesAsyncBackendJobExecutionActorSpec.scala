package cromwell.backend.impl.tes

import common.mock.MockSugar
import cromwell.core.logging.JobLogger
import cromwell.core.path.NioPath
import cromwell.filesystems.blob.{BlobFileSystemManager, BlobPath, WSMBlobSasTokenGenerator}
import org.mockito.ArgumentMatchers.any
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.time.Duration
import java.time.temporal.ChronoUnit
import scala.util.{Failure, Try}

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

  def generateMockWsmTokenGenerator: WSMBlobSasTokenGenerator = {
    val mockTokenGenerator = mock[WSMBlobSasTokenGenerator]
    val expectedTokenDuration: Duration = Duration.of(24, ChronoUnit.HOURS)
    mockTokenGenerator.getWSMSasFetchEndpoint(any[BlobPath], any[Option[Duration]]) returns Try(s"$testWsmEndpoint/api/workspaces/v1/$testWorkspaceId/resources/controlled/azure/storageContainer/$testContainerResourceId/getSasToken?sasExpirationDuration=${expectedTokenDuration.getSeconds.toInt}")
    mockTokenGenerator
  }
  def generateMockFsm: BlobFileSystemManager = {
    val mockFsm: BlobFileSystemManager = mock[BlobFileSystemManager]
    val mockGenerator: WSMBlobSasTokenGenerator = generateMockWsmTokenGenerator
    mockFsm.blobTokenGenerator returns mockGenerator
    mockFsm
  }
  //path to a blob file
  def generateMockBlobPath: BlobPath = {
    val mockBlobPath = mock[BlobPath]
    mockBlobPath.pathAsString returns someBlobUrl

    val mockFsm = generateMockFsm
    mockBlobPath.getFilesystemManager returns mockFsm

    val mockNioPath: NioPath = mock[NioPath]
    mockBlobPath.nioPath returns mockNioPath
    mockBlobPath
  }

  //Path to a file that isn't a blob file
  def generateMockDefaultPath: cromwell.core.path.Path = {
    val mockDefaultPath: cromwell.core.path.Path = mock[cromwell.core.path.Path]
    mockDefaultPath.pathAsString returns someNotBlobUrl
    mockDefaultPath
  }
  def pathGetter(pathString: String): Try[cromwell.core.path.Path] = {
    val mockBlob: BlobPath = generateMockBlobPath
    val mockDefault: cromwell.core.path.Path = generateMockDefaultPath
    if(pathString.contains(someBlobUrl)) Try(mockBlob) else Try(mockDefault)
  }

  def blobConverter(pathToConvert: Try[cromwell.core.path.Path]): Try[BlobPath] = {
    val mockBlob: BlobPath = generateMockBlobPath
    if(pathToConvert.get.pathAsString.contains(someBlobUrl)) Try(mockBlob) else Failure(new Exception("failed"))
  }

  it should "not return sas endpoint when no blob paths are provided" in {
    val mockLogger: JobLogger = mock[JobLogger]
    val emptyInputs: List[Input] = List()
    val bloblessInputs: List[Input] = List(notBlobInput_1, notBlobInput_2)
    TesAsyncBackendJobExecutionActor.determineWSMSasEndpointFromInputs(emptyInputs, pathGetter, mockLogger, blobConverter).isFailure shouldBe true
    TesAsyncBackendJobExecutionActor.determineWSMSasEndpointFromInputs(bloblessInputs, pathGetter, mockLogger, blobConverter).isFailure shouldBe true
  }

  it should "return a sas endpoint based on inputs when blob paths are provided" in {
    val mockLogger: JobLogger = mock[JobLogger]
    val expectedTokenLifetimeSeconds = 24 * 60 * 60 //assert that cromwell asks for 24h token duration.
    val expected = s"$testWsmEndpoint/api/workspaces/v1/$testWorkspaceId/resources/controlled/azure/storageContainer/$testContainerResourceId/getSasToken?sasExpirationDuration=${expectedTokenLifetimeSeconds}"
    val blobInput: List[Input] = List(blobInput_0)
    val blobInputs: List[Input] = List(blobInput_0, blobInput_1)
    val mixedInputs: List[Input] = List(notBlobInput_1, blobInput_0, blobInput_1)
    TesAsyncBackendJobExecutionActor.determineWSMSasEndpointFromInputs(blobInput, pathGetter, mockLogger, blobConverter).get shouldEqual expected
    TesAsyncBackendJobExecutionActor.determineWSMSasEndpointFromInputs(blobInputs, pathGetter, mockLogger, blobConverter).get shouldEqual expected
    TesAsyncBackendJobExecutionActor.determineWSMSasEndpointFromInputs(mixedInputs, pathGetter, mockLogger,  blobConverter).get shouldEqual expected
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
    val echoCommandSubstring = s"""echo "Saving sas token: $${$mockEnvironmentVariableNameFromWom:0:4}**** to environment variable $mockEnvironmentVariableNameFromWom...""""
    val generatedBashScript = TesAsyncBackendJobExecutionActor.generateLocalizedSasScriptPreamble(mockEnvironmentVariableNameFromWom, expectedEndpoint)

    generatedBashScript should include (beginSubstring)
    generatedBashScript should include (endSubstring)
    generatedBashScript should include (curlCommandSubstring)
    generatedBashScript should include (echoCommandSubstring)
    generatedBashScript should include (exportCommandSubstring)
  }
}
