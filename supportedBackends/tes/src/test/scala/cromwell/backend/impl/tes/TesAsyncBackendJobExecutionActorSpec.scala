package cromwell.backend.impl.tes

import common.mock.MockSugar
import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.standard.StandardAsyncJob
import cromwell.backend.{BackendJobDescriptorKey, BackendSpec}
import cromwell.core.TestKitSuite
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, NioPath}
import cromwell.filesystems.blob.{BlobFileSystemManager, BlobPath, WSMBlobSasTokenGenerator}
import cromwell.filesystems.http.HttpPathBuilder
import org.mockito.ArgumentMatchers.any
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import wom.graph.CommandCallNode

import java.time.Duration
import java.time.temporal.ChronoUnit
import java.util.UUID
import scala.concurrent.{ExecutionContext, Future}
import scala.util.{Failure, Try}

class TesAsyncBackendJobExecutionActorSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with BackendSpec
    with MockSugar
    with TableDrivenPropertyChecks {

  behavior of "TesAsyncBackendJobExecutionActor"

  type StandardAsyncRunInfo = Any
  type StandardAsyncPendingExecutionHandle =
    PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, TesRunStatus]
  implicit private val ec: ExecutionContext = system.dispatcher

  val fullyQualifiedName = "this.name.is.more.than.qualified"
  val workflowName = "mockWorkflow"
  val someBlobUrl =
    "https://lz813a3d637adefec2c6e88f.blob.core.windows.net/sc-d8143fd8-aa07-446d-9ba0-af72203f1794/nyxp6c/tes-internal/configuration/supported-vm-sizes"
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
  index = index + 1

  val blobInput_1 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someBlobUrl),
    path = someBlobUrl,
    `type` = Option("FILE"),
    content = None
  )
  index = index + 1

  val notBlobInput_1 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someNotBlobUrl + index),
    path = someNotBlobUrl + index,
    `type` = Option("FILE"),
    content = None
  )
  index = index + 1

  val notBlobInput_2 = Input(
    name = Option(fullyQualifiedName + "." + index),
    description = Option(workflowName + "." + fullyQualifiedName + "." + index),
    url = Option(someNotBlobUrl + index),
    path = someNotBlobUrl + index,
    `type` = Option("FILE"),
    content = None
  )

  val mockTesTaskLog = TaskLog(
    start_time = Option("2024-04-04T20:20:32.240066+00:00"),
    end_time = Option("2024-04-04T20:22:32.077818+00:00"),
    metadata = Option(Map("vm_price_per_hour_usd" -> "0.203")),
    logs = Option(mock[Seq[ExecutorLog]]),
    outputs = Option(mock[Seq[OutputFileLog]]),
    system_logs = Option(Seq("an error!"))
  )

  val mockTaskLog_1 = Task(
    Option(""),
    Option("Running"),
    Option("name"),
    Option("description"),
    Option(Seq.empty),
    Option(Seq.empty),
    null,
    null,
    null,
    null,
    Option(Seq(mockTesTaskLog))
  )

  val mockTaskLog_2 = Task(
    Option(""),
    Option("EXECUTOR_ERROR"),
    Option("name"),
    Option("description"),
    Option(Seq.empty),
    Option(Seq.empty),
    null,
    null,
    null,
    null,
    Option(Seq(mockTesTaskLog))
  )

  val mockMinimalView = MinimalTaskView("foo id", "Running")

  val mockMinimalView_2 = MinimalTaskView("abc123", "SYSTEM_ERROR")

  // Mock blob path functionality.
  val testWsmEndpoint = "https://wsm.mock.com/endpoint"
  val testWorkspaceId = "e58ed763-928c-4155-0000-fdbaaadc15f3"
  val testContainerResourceId = "e58ed763-928c-4155-1111-fdbaaadc15f3"

  def generateMockWsmTokenGenerator: WSMBlobSasTokenGenerator = {
    val mockTokenGenerator = mock[WSMBlobSasTokenGenerator]
    val expectedTokenDuration: Duration = Duration.of(24, ChronoUnit.HOURS)
    mockTokenGenerator.getWSMSasFetchEndpoint(any[BlobPath], any[Option[Duration]]) returns Try(
      s"$testWsmEndpoint/api/workspaces/v1/$testWorkspaceId/resources/controlled/azure/storageContainer/$testContainerResourceId/getSasToken?sasExpirationDuration=${expectedTokenDuration.getSeconds.toInt}"
    )
    mockTokenGenerator
  }
  def generateMockFsm: BlobFileSystemManager = {
    val mockFsm: BlobFileSystemManager = mock[BlobFileSystemManager]
    val mockGenerator: WSMBlobSasTokenGenerator = generateMockWsmTokenGenerator
    mockFsm.blobTokenGenerator returns mockGenerator
    mockFsm
  }
  // path to a blob file
  def generateMockBlobPath: BlobPath = {
    val mockBlobPath = mock[BlobPath]
    mockBlobPath.pathAsString returns someBlobUrl

    val mockFsm = generateMockFsm
    mockBlobPath.getFilesystemManager returns mockFsm

    val mockNioPath: NioPath = mock[NioPath]
    mockBlobPath.nioPath returns mockNioPath
    mockBlobPath
  }

  // Path to a file that isn't a blob file
  def generateMockDefaultPath: cromwell.core.path.Path = {
    val mockDefaultPath: cromwell.core.path.Path = mock[cromwell.core.path.Path]
    mockDefaultPath.pathAsString returns someNotBlobUrl
    mockDefaultPath
  }
  def pathGetter(pathString: String): Try[cromwell.core.path.Path] = {
    val mockBlob: BlobPath = generateMockBlobPath
    val mockDefault: cromwell.core.path.Path = generateMockDefaultPath
    if (pathString.contains(someBlobUrl)) Try(mockBlob) else Try(mockDefault)
  }

  def mockGetTaskLogs(): Future[Option[TaskLog]] =
    Future.successful(Option(mockTesTaskLog))

  def mockGetErrorLogs(handle: StandardAsyncPendingExecutionHandle): Future[Seq[String]] = {
    val logs = mockTesTaskLog
    val systemLogs = logs.system_logs.get
    Future.successful(systemLogs)
  }

  def mockFetchFullTaskView_1(handle: StandardAsyncPendingExecutionHandle): Future[Task] =
    Future.successful(mockTaskLog_1)

  def mockFetchFullTaskView_2(handle: StandardAsyncPendingExecutionHandle): Future[Task] =
    Future.successful(mockTaskLog_2)

  def mockFetchMinimalTaskView(handle: StandardAsyncPendingExecutionHandle): Future[MinimalTaskView] =
    Future.successful(mockMinimalView)

  def mockFetchMinimalTaskView_2(handle: StandardAsyncPendingExecutionHandle): Future[MinimalTaskView] =
    Future.successful(mockMinimalView_2)

  def mockGetTesStatus(state: Option[String], withCostData: Option[TesVmCostData], jobId: String): TesRunStatus =
    state match {
      case s if s.contains("COMPLETE") =>
        Complete(withCostData)

      case s if s.contains("CANCELED") =>
        Cancelled(withCostData)

      case s if s.contains("EXECUTOR_ERROR") =>
        Failed(Seq.empty[String], withCostData)

      case s if s.contains("SYSTEM_ERROR") =>
        Error()

      case _ =>
        Running(withCostData)
    }

  def mockTellMetadata(metadataMap: Map[String, Any]): Unit = {
    Future.successful(metadataMap)
    ()
  }

  def blobConverter(pathToConvert: Try[cromwell.core.path.Path]): Try[BlobPath] = {
    val mockBlob: BlobPath = generateMockBlobPath
    if (pathToConvert.get.pathAsString.contains(someBlobUrl)) Try(mockBlob) else Failure(new Exception("failed"))
  }

  it should "not return sas endpoint when no blob paths are provided" in {
    val mockLogger: JobLogger = mock[JobLogger]
    val emptyInputs: List[Input] = List()
    val bloblessInputs: List[Input] = List(notBlobInput_1, notBlobInput_2)
    TesAsyncBackendJobExecutionActor
      .determineWSMSasEndpointFromInputs(emptyInputs, pathGetter, mockLogger, blobConverter)
      .isFailure shouldBe true
    TesAsyncBackendJobExecutionActor
      .determineWSMSasEndpointFromInputs(bloblessInputs, pathGetter, mockLogger, blobConverter)
      .isFailure shouldBe true
  }

  it should "return a sas endpoint based on inputs when blob paths are provided" in {
    val mockLogger: JobLogger = mock[JobLogger]
    val expectedTokenLifetimeSeconds = 24 * 60 * 60 // assert that cromwell asks for 24h token duration.
    val expected =
      s"$testWsmEndpoint/api/workspaces/v1/$testWorkspaceId/resources/controlled/azure/storageContainer/$testContainerResourceId/getSasToken?sasExpirationDuration=${expectedTokenLifetimeSeconds}"
    val blobInput: List[Input] = List(blobInput_0)
    val blobInputs: List[Input] = List(blobInput_0, blobInput_1)
    val mixedInputs: List[Input] = List(notBlobInput_1, blobInput_0, blobInput_1)
    TesAsyncBackendJobExecutionActor
      .determineWSMSasEndpointFromInputs(blobInput, pathGetter, mockLogger, blobConverter)
      .get shouldEqual expected
    TesAsyncBackendJobExecutionActor
      .determineWSMSasEndpointFromInputs(blobInputs, pathGetter, mockLogger, blobConverter)
      .get shouldEqual expected
    TesAsyncBackendJobExecutionActor
      .determineWSMSasEndpointFromInputs(mixedInputs, pathGetter, mockLogger, blobConverter)
      .get shouldEqual expected
  }

  it should "contain expected strings in the bash script" in {
    val mockEnvironmentVariableNameFromWom = "mock_env_var_for_storing_sas_token"
    val expectedEndpoint =
      s"$testWsmEndpoint/api/workspaces/v1/$testWorkspaceId/resources/controlled/azure/storageContainer/$testContainerResourceId/getSasToken"

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
         |                    -H "Authorization: Bearer $${BEARER_TOKEN}" \\
         |                    -H "Content-Length: 0" \\
         |                    -d "")
         |""".stripMargin
    val exportCommandSubstring =
      s"""export $mockEnvironmentVariableNameFromWom=$$(echo "$${sas_response_json}" | jq -r '.token')"""
    val echoCommandSubstring =
      s"""echo "Saving sas token: $${$mockEnvironmentVariableNameFromWom:0:4}**** to environment variable $mockEnvironmentVariableNameFromWom...""""
    val generatedBashScript =
      TesAsyncBackendJobExecutionActor.generateLocalizedSasScriptPreamble(mockEnvironmentVariableNameFromWom,
                                                                          expectedEndpoint
      )

    generatedBashScript should include(beginSubstring)
    generatedBashScript should include(endSubstring)
    generatedBashScript should include(curlCommandSubstring)
    generatedBashScript should include(echoCommandSubstring)
    generatedBashScript should include(exportCommandSubstring)
  }

  it should "return expected task end time" in {
    val endTime = TesAsyncBackendJobExecutionActor.getTaskEndTime(mockGetTaskLogs())
    whenReady(endTime) { m =>
      m.get shouldBe "2024-04-04T20:22:32.077818+00:00"
    }
  }

  it should "return error logs for error state" in {
    val errorLogs = TesAsyncBackendJobExecutionActor.getErrorSeq(mockGetTaskLogs())

    whenReady(errorLogs) { m =>
      m.get shouldBe Seq("an error!")
    }
  }

  it should "return expected status determined by cost data" in {
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new StandardAsyncPendingExecutionHandle(null, runId, None, None)
    val tesStatusWithData = TesAsyncBackendJobExecutionActor.queryStatusAndMaybeCostData(handle,
                                                                                         true,
                                                                                         mockFetchFullTaskView_1,
                                                                                         mockFetchMinimalTaskView,
                                                                                         mockGetTesStatus,
                                                                                         mockTellMetadata
    )

    whenReady(tesStatusWithData) { s =>
      s shouldEqual (Running(
        Option(
          TesVmCostData(Option("2024-04-04T20:20:32.240066+00:00"), None, Option("0.203"))
        )
      ))
    }

    val tesStatusNoData = TesAsyncBackendJobExecutionActor.queryStatusAndMaybeCostData(handle,
                                                                                       false,
                                                                                       mockFetchFullTaskView_1,
                                                                                       mockFetchMinimalTaskView,
                                                                                       mockGetTesStatus,
                                                                                       mockTellMetadata
    )

    whenReady(tesStatusNoData) { s =>
      s shouldEqual (Running(None))
    }
  }

  it should "return expected error states with expected data" in {
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val costData = Option(TesVmCostData(Option("time"), Option("time"), Option("0.203")))
    val handle = new StandardAsyncPendingExecutionHandle(null, runId, None, Option(Running(None)))

    // No cost data
    val actualStatus = TesAsyncBackendJobExecutionActor.pollTesStatus(handle,
                                                                      true,
                                                                      mockFetchFullTaskView_2,
                                                                      mockFetchMinimalTaskView,
                                                                      mockGetTesStatus,
                                                                      mockTellMetadata,
                                                                      mockGetErrorLogs
    )
    val expectedStatus = Failed(List("an error!"), None)

    whenReady(actualStatus) { s =>
      s shouldBe expectedStatus
    }

    // With cost data on handle
    val handleWithData = new StandardAsyncPendingExecutionHandle(null, runId, None, Option(Running(costData)))

    val actualStatusWithData = TesAsyncBackendJobExecutionActor.pollTesStatus(handleWithData,
                                                                              true,
                                                                              mockFetchFullTaskView_2,
                                                                              mockFetchMinimalTaskView,
                                                                              mockGetTesStatus,
                                                                              mockTellMetadata,
                                                                              mockGetErrorLogs
    )
    val expectedStatus2 = Failed(List("an error!"), costData)

    whenReady(actualStatusWithData) { s =>
      s shouldBe expectedStatus2
    }
  }

  it should "return expected error states with expected data false" in {
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val costData = Option(TesVmCostData(Option("time"), Option("time"), Option("0.203")))
    val handle = new StandardAsyncPendingExecutionHandle(null, runId, None, Option(Running(costData)))
    val actualStatus = TesAsyncBackendJobExecutionActor.pollTesStatus(handle,
                                                                      false,
                                                                      mockFetchFullTaskView_2,
                                                                      mockFetchMinimalTaskView_2,
                                                                      mockGetTesStatus,
                                                                      mockTellMetadata,
                                                                      mockGetErrorLogs
    )
    val expectedStatus = Error(List("an error!"), costData)

    whenReady(actualStatus) { s =>
      s shouldBe expectedStatus
    }
  }

  private val httpPathTestCases = Table(
    ("test name", "http path", "local path in input dir"),
    (
      "strip simple kv query params",
      "http://example.com/my_sample.bam?k1=v1&k2=v2",
      "example.com/my_sample.bam"
    ),
    (
      "handle http paths without query params",
      "http://example.com/my_sample.bam",
      "example.com/my_sample.bam"
    ),
    (
      "handle http paths without params but with a ?",
      "http://example.com/my_sample.bam?",
      "example.com/my_sample.bam"
    ),
    (
      "handle a blob file with SAS token attached",
      "https://lzbc096764ae93ffff9f406e.blob.core.windows.net/sc-a7f7a9e0-2dcf-465c-997b-a276090a52da/workspace-services/cbas/terra-app-2f577477-763b-4e27-8e28-b03d91b6f3be/cromwell-workflow-logs/workflow.c621a5df-37f1-422d-b91a-1a65f6112a6a.log?sv=2023-11-03&spr=https&st=2024-04-09T23%3A35%3A37Z&se=2024-04-10T07%3A50%3A37Z&sr=c&sp=racwdlt&sig=REDACTEDS&rscd=100067995116984528334",
      "lzbc096764ae93ffff9f406e.blob.core.windows.net/sc-a7f7a9e0-2dcf-465c-997b-a276090a52da/workspace-services/cbas/terra-app-2f577477-763b-4e27-8e28-b03d91b6f3be/cromwell-workflow-logs/workflow.c621a5df-37f1-422d-b91a-1a65f6112a6a.log"
    ),
    (
      "handle an http path with fragment",
      "http://example.com/my_sample.bam#my_favorite_part",
      "example.com/my_sample.bam"
    ),
    (
      "handle an http path with fragment and query params",
      "http://example.com/my_sample.bam?k=yourface#my_favorite_part",
      "example.com/my_sample.bam"
    )
  )

  forAll(httpPathTestCases) { (testName, httpPath, localPathInInputDir) =>
    it should testName in {
      val wd = buildWdlWorkflowDescriptor(TestWorkflows.HelloWorld)
      val call: CommandCallNode = wd.callable.taskCallNodes.head
      val jobKey = BackendJobDescriptorKey(call, None, 1)
      val jobPaths = TesJobPaths(jobKey, wd, TesTestConfig.backendConfig)
      val commandDirectory = DefaultPathBuilder.build("/my/command/dir").get
      val httpBuilder = new HttpPathBuilder()

      val httpPathWithParams = httpBuilder.build(httpPath)
      val actual = TesAsyncBackendJobExecutionActor.mapInputPath(httpPathWithParams.get, jobPaths, commandDirectory)
      actual shouldBe s"${jobPaths.callInputsDockerRoot}/$localPathInInputDir"
    }
  }
}
