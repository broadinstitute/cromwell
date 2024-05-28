package cromwell.backend.impl.tes

import akka.actor.ActorRef
import akka.http.scaladsl.client.RequestBuilding.Get
import akka.http.scaladsl.model.{ContentTypes, HttpEntity, HttpProtocols, HttpResponse, StatusCodes}
import akka.testkit.ImplicitSender
import com.typesafe.config.Config
import common.mock.MockSugar
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.{
  BackendJobDescriptor,
  BackendJobDescriptorKey,
  BackendSpec,
  MinimumRuntimeSettings,
}
import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.io.JobPathsSpecHelper.DummyStandardPaths
import cromwell.backend.standard.{
  DefaultStandardAsyncExecutionActorParams,
  StandardAsyncExecutionActorParams,
  StandardAsyncJob,
  StandardExpressionFunctions,
  StandardExpressionFunctionsParams
}
import cromwell.core.{CallContext, TestKitSuite}
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, NioPath, PathBuilder}
import cromwell.filesystems.blob.{BlobFileSystemManager, BlobPath, WSMBlobSasTokenGenerator}
import cromwell.filesystems.http.HttpPathBuilder
import org.mockito.ArgumentMatchers.any
import org.mockito.Mockito.when
import org.scalatest.{BeforeAndAfter, PrivateMethodTester}
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import org.slf4j.Logger
import spray.json.DefaultJsonProtocol
import wom.graph.CommandCallNode

import java.net.URI
import java.time.Duration
import java.time.temporal.ChronoUnit
import java.util.UUID
import scala.concurrent.{ExecutionContext, ExecutionContextExecutor, Future, Promise}
import scala.util.{Failure, Try}

class TesAsyncBackendJobExecutionActorSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with BackendSpec
    with ImplicitSender
    with BeforeAndAfter
    with MockSugar
    with DefaultJsonProtocol
    with PrivateMethodTester
    with TableDrivenPropertyChecks {

  val YoSup: String =
    s"""
       |task sup {
       |  String addressee
       |  command {
       |    echo "yo sup $${addressee}!"
       |  }
       |  output {
       |    String salutation = read_string(stdout())
       |  }
       |  runtime {
       |    docker: "alpine:latest"
       |    queueArn: "arn:aws:batch:us-east-1:111222333444:job-queue/job-queue"
       |  }
       |}
       |
       |workflow wf_sup {
       |  call sup
       |}
    """.stripMargin

  lazy val mockPathBuilderLocal: PathBuilder = DefaultPathBuilder

  private lazy val TestableCallContext =
    CallContext(mockPathBuilderLocal.build("/root").get, DummyStandardPaths, isDocker = false)

  lazy val TestableStandardExpressionFunctionsParams: StandardExpressionFunctionsParams =
    new StandardExpressionFunctionsParams {
      override lazy val pathBuilders: List[PathBuilder] = List(mockPathBuilderLocal)
      override lazy val callContext: CallContext = TestableCallContext
      override val ioActorProxy: ActorRef = simpleIoActor
      override val executionContext: ExecutionContext = system.dispatcher
    }
  lazy val TestableTesExpressionFunctions: TesExpressionFunctions = new TesExpressionFunctions(
    TestableStandardExpressionFunctionsParams
  )
  private def buildInitializationData(jobDescriptor: BackendJobDescriptor,
                                      tesConfiguration: TesConfiguration,
                                      config: Config
  ) = {
    val workflowPaths = TesWorkflowPaths(
      jobDescriptor.workflowDescriptor,
      config
    )
    val runtimeAttributesBuilder = TesRuntimeAttributes.runtimeAttributesBuilder(Option(config))
    TesBackendInitializationData(workflowPaths, runtimeAttributesBuilder, tesConfiguration)
  }

  class TestableTesAsyncBackendJobExecutionActor(params: StandardAsyncExecutionActorParams,
                                                 functions: StandardExpressionFunctions
  ) extends TesAsyncBackendJobExecutionActor(params) {

    type StandardAsyncRunInfo
    type StandardAsyncRunState
    type TesPendingExecutionHandle =
      PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState]

    def this(jobDescriptor: BackendJobDescriptor,
             promise: Promise[BackendJobExecutionResponse],
             configuration: TesConfiguration,
             config: Config,
             functions: StandardExpressionFunctions = TestableTesExpressionFunctions,
             singletonActor: ActorRef = emptyActor,
             ioActor: ActorRef = mockIoActor
    ) =
      this(
        DefaultStandardAsyncExecutionActorParams(
          jobIdKey = TesAsyncBackendJobExecutionActor.JobIdKey,
          serviceRegistryActor = singletonActor,
          ioActor = ioActor,
          jobDescriptor = jobDescriptor,
          configurationDescriptor = configuration.configurationDescriptor,
          backendInitializationDataOption = Option(buildInitializationData(jobDescriptor, configuration, config)),
          backendSingletonActorOption = Option(singletonActor),
          completionPromise = promise,
          minimumRuntimeSettings = MinimumRuntimeSettings()
        ),
        functions
      )

    override lazy val jobLogger: JobLogger = new JobLogger(
      "TestLogger",
      workflowIdForLogging,
      rootWorkflowIdForLogging,
      jobTag,
      akkaLogger = Option(log)
    ) {
      override def tag: String = s"$name [UUID(${workflowIdForLogging.shortString})$jobTag]"
      override val slf4jLoggers: Set[Logger] = Set.empty
    }

    override lazy val backendEngineFunctions: StandardExpressionFunctions = functions
  }

  behavior of "TesAsyncBackendJobExecutionActor"

  // private val timeout = 60 seconds
  type StandardAsyncRunInfo = Any
  type StandardAsyncRunState = TesRunStatus
  type StandardAsyncPendingExecutionHandle =
    PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState]
  implicit private val ec: ExecutionContextExecutor = system.dispatcher
  val standardParams: StandardAsyncExecutionActorParams = mock[StandardAsyncExecutionActorParams]
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

  val mockTaskLog = Task(
    Option(""), Option("SYSTEM_ERROR"), Option("name"), Option("description"), Option(Seq.empty), Option(Seq.empty), null, null, null, null, Option(Seq(mockTesTaskLog))
  )

  val mockMinimalView = MinimalTaskView("foo id", "Running")

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

  def mockGetTaskLogs(): Future[TaskLog] =
    Future.successful(mockTesTaskLog)

  def mockGetErrorLogs(handle: StandardAsyncPendingExecutionHandle, handler: HttpHandler): Future[Seq[String]] = {
    val logs = mockTesTaskLog
    val systemLogs = logs.system_logs.get
    Future.successful(systemLogs)
  }

  def mockQueryAndCostData(handle: StandardAsyncPendingExecutionHandle,
                           data: Boolean,
                           handler: HttpHandler
  ): Future[TesRunStatus] = {
    val fetchedCostData = Option(TesVmCostData(Option(""), Option("0,203")))
    Future.successful(Running(fetchedCostData))
  }

  def mockFetchFullTaskView(handle: StandardAsyncPendingExecutionHandle, handler: HttpHandler): Future[Task] =
    Future.successful(mockTaskLog)

  def mockFetchMinimalTaskView(handle: StandardAsyncPendingExecutionHandle,
                               handler: HttpHandler
  ): Future[MinimalTaskView] =
    Future.successful(mockMinimalView)

  def mockGetTesStatus(state: Option[String], withCostData: Option[TesVmCostData], jobId: String): TesRunStatus = {
    state match {
      case s if s.contains("COMPLETE") =>
        Complete(withCostData)

      case s if s.contains("CANCELED") =>
        Cancelled()

      case s if s.contains("EXECUTOR_ERROR") =>
        Failed()

      case s if s.contains("SYSTEM_ERROR") =>
        System.out.print("TES STATUS" + s)
        Error()

      case _ => Running(withCostData)
    }
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
    val mockHttpClient = mock[HttpHandler]
    val mockUri = new URI("http://example.com")
    val expectedGet = Get(mockUri.toString)
    val mockTesResponse = new HttpResponse(
      StatusCodes.OK,
      null,
      HttpEntity(ContentTypes.`application/json`, tesResponse),
      HttpProtocols.`HTTP/1.1`
    )
    val httpSuccess = Future.successful(mockTesResponse)

    val hihi = mockHttpClient.singleRequest(expectedGet, null)
    when(hihi).thenReturn(httpSuccess)

    val tesStatusWithData = TesAsyncBackendJobExecutionActor.queryStatusAndMaybeCostData(handle,
      true,
                                                                                 mockHttpClient,
                                                                                 mockFetchFullTaskView,
                                                                                 mockFetchMinimalTaskView,
                                                                                 mockGetTesStatus,
                                                                                 mockTellMetadata
    )

    whenReady(tesStatusWithData) { s =>
      s shouldEqual(Running(Option(TesVmCostData(Option("2024-04-04T20:20:32.240066+00:00"), Option("0.203")))))
    }

    val tesStatusNoData = TesAsyncBackendJobExecutionActor.queryStatusAndMaybeCostData(handle,
      false,
      mockHttpClient,
      mockFetchFullTaskView,
      mockFetchMinimalTaskView,
      mockGetTesStatus,
      mockTellMetadata
    )

    whenReady(tesStatusNoData) { s =>
      s shouldEqual(Running(None))
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

  val tesResponse =
    """{
      |  "id": "0166e1e8_2a12c2e2e86a44cd88c6c41322f85c9a",
      |  "state": "RUNNING",
      |  "name": "fetch_sra_to_bam.Fetch_SRA_to_BAM",
      |  "description": "0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6:BackendJobDescriptorKey_CommandCallNode_fetch_sra_to_bam.Fetch_SRA_to_BAM:-1:1",
      |  "inputs": [
      |    {
      |      "name": "commandScript",
      |      "description": "fetch_sra_to_bam.Fetch_SRA_to_BAM.commandScript",
      |      "url": "https://lz813a3d637adefec2c6e88f.blob.core.windows.net/sc-5fa76398-ac93-4d16-944f-15e844f79e7b/workspace-services/cbas/terra-app-2aece314-9a52-4b70-86d8-5d6f7d3f2189/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/script",
      |      "path": "/cromwell-executions/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/script",
      |      "type": "FILE",
      |      "content": null,
      |      "streamable": false
      |    }
      |  ],
      |  "outputs": [
      |    {
      |      "name": "rc",
      |      "description": "fetch_sra_to_bam.Fetch_SRA_to_BAM.rc",
      |      "url": "https://lz813a3d637adefec2c6e88f.blob.core.windows.net/sc-5fa76398-ac93-4d16-944f-15e844f79e7b/workspace-services/cbas/terra-app-2aece314-9a52-4b70-86d8-5d6f7d3f2189/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/rc",
      |      "path": "/cromwell-executions/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/rc",
      |      "type": "FILE"
      |    },
      |    {
      |      "name": "stdout",
      |      "description": "fetch_sra_to_bam.Fetch_SRA_to_BAM.stdout",
      |      "url": "https://lz813a3d637adefec2c6e88f.blob.core.windows.net/sc-5fa76398-ac93-4d16-944f-15e844f79e7b/workspace-services/cbas/terra-app-2aece314-9a52-4b70-86d8-5d6f7d3f2189/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/stdout",
      |      "path": "/cromwell-executions/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/stdout",
      |      "type": "FILE"
      |    },
      |    {
      |      "name": "stderr",
      |      "description": "fetch_sra_to_bam.Fetch_SRA_to_BAM.stderr",
      |      "url": "https://lz813a3d637adefec2c6e88f.blob.core.windows.net/sc-5fa76398-ac93-4d16-944f-15e844f79e7b/workspace-services/cbas/terra-app-2aece314-9a52-4b70-86d8-5d6f7d3f2189/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/stderr",
      |      "path": "/cromwell-executions/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/stderr",
      |      "type": "FILE"
      |    },
      |    {
      |      "name": "commandScript",
      |      "description": "fetch_sra_to_bam.Fetch_SRA_to_BAM.commandScript",
      |      "url": "https://lz813a3d637adefec2c6e88f.blob.core.windows.net/sc-5fa76398-ac93-4d16-944f-15e844f79e7b/workspace-services/cbas/terra-app-2aece314-9a52-4b70-86d8-5d6f7d3f2189/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/script",
      |      "path": "/cromwell-executions/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/script",
      |      "type": "FILE"
      |    },
      |    {
      |      "name": "fetch_sra_to_bam.Fetch_SRA_to_BAM.output.0",
      |      "description": "fetch_sra_to_bam.Fetch_SRA_to_BAM.output.0",
      |      "url": "https://lz813a3d637adefec2c6e88f.blob.core.windows.net/sc-5fa76398-ac93-4d16-944f-15e844f79e7b/workspace-services/cbas/terra-app-2aece314-9a52-4b70-86d8-5d6f7d3f2189/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/OUT_PLATFORM",
      |      "path": "/cromwell-executions/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/OUT_PLATFORM",
      |      "type": "FILE"
      |    },
      |    {
      |      "name": "fetch_sra_to_bam.Fetch_SRA_to_BAM.output.1",
      |      "description": "fetch_sra_to_bam.Fetch_SRA_to_BAM.output.1",
      |      "url": "https://lz813a3d637adefec2c6e88f.blob.core.windows.net/sc-5fa76398-ac93-4d16-944f-15e844f79e7b/workspace-services/cbas/terra-app-2aece314-9a52-4b70-86d8-5d6f7d3f2189/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/OUT_STRAIN",
      |      "path": "/cromwell-executions/fetch_sra_to_bam/0166e1e8-59c6-4ff7-abe9-4aa2e20f24f6/call-Fetch_SRA_to_BAM/execution/OUT_STRAIN",
      |      "type": "FILE"
      |    },
      |  ],
      |  "logs": [
      |    {
      |      "logs": [
      |        {
      |          "start_time": "2024-04-04T20:20:32.240066+00:00",
      |          "end_time": "2024-04-04T20:22:32.077818+00:00",
      |          "stdout": null,
      |          "stderr": null,
      |          "exit_code": 0
      |        }
      |      ],
      |      "metadata": {
      |        "vm_size": "Standard_D5_v2",
      |        "vm_family": "standardDv2Family",
      |        "vm_low_priority": "true",
      |        "vm_memory_in_gib": "56.0",
      |        "vm_vcpus_available": "16",
      |        "vm_price_per_hour_usd": "0.203",
      |        "vm_hyper_v_generations": "[\"V1\"]",
      |        "vm_max_data_disk_count": "64",
      |        "encryption_at_host_supported": "false",
      |        "vm_resource_disk_size_in_gib": "800.0",
      |        "cromwell_rc": "0"
      |      },
      |      "start_time": "2024-04-04T20:17:01.0408443+00:00",
      |      "end_time": "2024-04-04T20:22:36.6841675+00:00",
      |      "outputs": [],
      |      "system_logs": []
      |    }
      |  ],
      |  "creation_time": "2024-04-04T20:16:50.3092382+00:00"
      |}""".stripMargin
}
