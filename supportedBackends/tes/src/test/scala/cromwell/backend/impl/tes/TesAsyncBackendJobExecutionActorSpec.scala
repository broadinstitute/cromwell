package cromwell.backend.impl.tes

import akka.actor.{ActorRef, Props}
import akka.http.scaladsl.client.RequestBuilding.Get
import akka.testkit.{ImplicitSender, TestActorRef}
import com.typesafe.config.Config
import common.mock.MockSugar
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.{BackendJobDescriptor, BackendJobDescriptorKey, BackendSpec, BackendWorkflowDescriptor, MinimumRuntimeSettings, RuntimeAttributeDefinition}
import cromwell.backend.async.PendingExecutionHandle
import cromwell.backend.impl.tes.TesTestConfig.TesBackendConfigurationDescriptor
import cromwell.backend.io.JobPathsSpecHelper.DummyStandardPaths
import cromwell.backend.standard.{DefaultStandardAsyncExecutionActorParams, StandardAsyncExecutionActorParams, StandardAsyncJob, StandardExpressionFunctions, StandardExpressionFunctionsParams}
import cromwell.core.callcaching.NoDocker
import cromwell.core.labels.Labels
import cromwell.core.{CallContext, HogGroup, TestKitSuite, WorkflowId, WorkflowOptions}
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, NioPath, PathBuilder}
import cromwell.filesystems.blob.{BlobFileSystemManager, BlobPath, WSMBlobSasTokenGenerator}
import cromwell.filesystems.http.HttpPathBuilder
import cromwell.services.keyvalue.KeyValueServiceActor.KvResponse
import cromwell.util.JsonFormatting.WomValueJsonFormatter.WomValueJsonFormat
import org.mockito.ArgumentMatchers.any
import org.scalatest.{BeforeAndAfter, PrivateMethodTester}
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks
import org.slf4j.Logger
import spray.json.{DefaultJsonProtocol, JsObject, JsValue, enrichAny}
import wdl.draft2.model.{Draft2ImportResolver, FullyQualifiedName, WdlNamespaceWithWorkflow}
import wom.expression.NoIoFunctionSet
import wom.graph.CommandCallNode
import wom.values.{WomString, WomValue}
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomExecutableMakers._
import wdl.transforms.draft2.wdlom2wom._
import wom.transforms.WomExecutableMaker.ops._
import wom.transforms.WomWorkflowDefinitionMaker.ops._

import java.time.Duration
import java.time.temporal.ChronoUnit
// import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, ExecutionContextExecutor, Future, Promise}
// import scala.language.postfixOps
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

  private val Inputs: Map[FullyQualifiedName, WomValue] = Map("wf_sup.sup.addressee" -> WomString("dog"))

  private val NoOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))

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
    type TesPendingExecutionHandle = PendingExecutionHandle[StandardAsyncJob, StandardAsyncRunInfo, StandardAsyncRunState]

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

  private val configuration = new TesConfiguration(TesBackendConfigurationDescriptor)
  private val config = mock[Config]
  private val runtimeAttributesBuilder = TesRuntimeAttributes.runtimeAttributesBuilder(Option(config))


  private def buildJobDescriptor(): BackendJobDescriptor = {
    val attempt = 1
    val wdlNamespace = WdlNamespaceWithWorkflow.load(YoSup, Seq.empty[Draft2ImportResolver]).get
    val womDefinition = wdlNamespace.workflow.toWomWorkflowDefinition(isASubworkflow = false)
      .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))

    wdlNamespace.toWomExecutable(Option(Inputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true) match {
      case Right(womExecutable) =>
        val inputs = for {
          combined <- womExecutable.resolvedExecutableInputs
          (port, resolvedInput) = combined
          value <- resolvedInput.select[WomValue]
        } yield port -> value

        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womDefinition,
          inputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val job = workflowDescriptor.callable.taskCallNodes.head
        val key = BackendJobDescriptorKey(job, None, attempt)
        val runtimeAttributes = makeRuntimeAttributes(job)
        val prefetchedKvEntries = Map[String, KvResponse]()
        BackendJobDescriptor(workflowDescriptor,
          key,
          runtimeAttributes,
          fqnWdlMapToDeclarationMap(Inputs),
          NoDocker,
          None,
          prefetchedKvEntries
        )
      case Left(badtimes) => fail(badtimes.toList.mkString(", "))
    }
  }

  def buildTestActorRef(functions: StandardExpressionFunctions): TestActorRef[TestableTesAsyncBackendJobExecutionActor] = {
    // For this test we say that all previous attempts were preempted:
    val jobDescriptor = buildJobDescriptor()
    val mockConfig = mock[Config]
    val props = Props(
      new TestableTesAsyncBackendJobExecutionActor(
        jobDescriptor,
        Promise(),
        configuration,
        mockConfig,
        functions,
        emptyActor,
        mockIoActor
      )
    )
    TestActorRef(props, s"TestableTesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")
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
    null,
    null,
    null,
    null,
    null,
    null,
    null,
    null,
    null,
    null,
    null,
  )

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

  def mockGetTaskLogs(): Future[TaskLog] = {
    Future.successful(mockTesTaskLog)
  }

  def mockGetErrorLogs(handle: StandardAsyncPendingExecutionHandle): Future[Seq[String]] = {
    val logs = mockTesTaskLog
    val systemLogs = logs.system_logs.get
    Future.successful(systemLogs)
  }

  def mockQueryAndCostData(handle: StandardAsyncPendingExecutionHandle, data: Boolean): Future[TesRunStatus] = {
    val fetchedCostData = Option(TesVmCostData(Option(""), Option("0,203")))
    if (data) {
      handle.previousState match {
        case Running(Some(_), Some(_)) => Future.successful(Running(fetchedCostData))
        case Running(None) => Future.successful(Running(fetchedCostData))
        case None => Future.successful(Running(fetchedCostData))
      }
    } else {
      handle.previousState match {
        case Running(Some(_), Some(_)) => Future.successful(Running(fetchedCostData))
        case Running(None, Some(_)) | Running(Some(_), None) =>
        case Running(None) => Future.successful(Running(None))
        case None => Future.successful(Running(None))
      }
    }
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

  it should "return and change my name" in {
    val mockHttpClient = mock[HttpHandler]
    val expectedGet = Get()
    val hihi = mockHttpClient.singleRequest(expectedGet, )
  }

//  it should "do something" in {
//    val actorRef = buildTestActorRef(TestableTesExpressionFunctions)
//    val runId = StandardAsyncJob(UUID.randomUUID().toString)
//    val handle = new TesPendingExecutionHandle(null, runId, None, None)// mock[StandardAsyncPendingExecutionHandle]
//    val mockTesVmCostData = TesVmCostData(Option("time"), Option("0.203"))
//    val status = Running(Option(mockTesVmCostData))
////    val expectedRequest = Get(mockUri.toString)
////    val mockHttpClient = mock[HttpClientBuilder]
//
////    when(mockHttpClient.build().execute(expectedRequest)).thenReturn(mockTesTaskLog)
//    val actor = actorRef.underlyingActor
//    // val yay = actorRef ! System.out.print("hello")
//     val result = actor.onTaskComplete(status, handle)
//    //functionshere.pollStatus(handle)
//    // when(functionshere.pollStatus(handle)).thenReturn(status)
//
//    // val result = Await.result(functionshere.pollStatusAsync(handle), timeout)
//    System.out.print(result)
//    //System.out.print(actorRef.underlyingActor)
//
//  }

//  it should "do something else" in {
//    val runId = StandardAsyncJob(UUID.randomUUID().toString)
//    val handle = new TesPendingExecutionHandle(null, runId, None, None)
////    val actor = TestActorRef[TestTesAsyncBackendJobExecutionActor]//buildTestActorRef(TestableTesExpressionFunctions)
//    val tesExecutionActor = system.actorOf(Props(new TestTesAsyncBackendJobExecutionActor(mock[StandardAsyncExecutionActorParams]))) //new TestTesAsyncBackendJobExecutionActor(mock[StandardAsyncExecutionActorParams])
//
//    val tesStatus = tesExecutionActor.pollTesStatus(handle, false)
//
//   val mockTesVmCostData = TesVmCostData(Option("time"), Option("0.203"))
//   val status = Running(Option(mockTesVmCostData))
//
//    // val pollTesStatus = PrivateMethod[Future[TesRunStatus]](Symbol("pollTesStatus"))
//
//        // testing a private method see https://www.scalatest.org/user_guide/using_PrivateMethodTester
//    // val pollStat = idk invokePrivate pollTesStatus(handle, false)
//    //when(mocker.pollStatusAsync(handle)).thenReturn(Future.successful(Running(Option(mockTesVmCostData))))
//    // when(pollStat).thenReturn(Future.successful(Running(Option(mockTesVmCostData))))
//  System.out.print(tesStatus)
//    // val webReq =
//    // when(actor.underlyingActor.fetchFullTaskView(handle)).thenReturn(Future.successful(mockTaskLog))
//    val mockmock = mock[TestableTesAsyncBackendJobExecutionActor]
//    doReturn(Future.successful(status)).when(mockmock).queryStatusAndMaybeCostData(handle, false)
////    doReturn(Future.fromTry(Success(mockTaskLog))).when(mockmock).fetchFullTaskView(handle)
//
////    when(actor.underlyingActor.queryStatusAndMaybeCostData(handle, false)).thenReturn(Future.successful(status))
//    //val query = actor.underlyingActor.queryStatusAndMaybeCostData(handle, false)
//    // System.out.print(query)
////    val result = actor.underlyingActor.pollTesStatus(handle, false)
////    System.out.print(result)
//
////    whenReady(query) { m =>
////      System.out.print(m)
//  }

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

  private def makeRuntimeAttributes(job: CommandCallNode) = {
    val evaluatedAttributes =
      RuntimeAttributeDefinition.evaluateRuntimeAttributes(job.callable.runtimeAttributes, null, Map.empty)
    RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributesBuilder.definitions.toSet, NoOptions)(
      evaluatedAttributes.getOrElse(fail("Failed to evaluate runtime attributes"))
    )
  }
}
