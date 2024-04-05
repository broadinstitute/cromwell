package cromwell.backend.google.pipelines.common

import java.nio.file.Paths
import java.util.UUID

import _root_.io.grpc.Status
import _root_.wdl.draft2.model._
import akka.actor.{ActorRef, Props}
import akka.testkit.{ImplicitSender, TestActorRef, TestDuration, TestProbe}
import cats.data.NonEmptyList
import com.google.api.client.http.HttpRequest
import com.google.cloud.NoCredentials
import common.collections.EnhancedCollections._
import cromwell.backend.BackendJobExecutionActor.{
  BackendJobExecutionResponse,
  JobFailedNonRetryableResponse,
  JobFailedRetryableResponse
}
import cromwell.backend._
import cromwell.backend.async.AsyncBackendJobExecutionActor.{Execute, ExecutionMode}
import cromwell.backend.async.{
  AbortedExecutionHandle,
  ExecutionHandle,
  FailedNonRetryableExecutionHandle,
  FailedRetryableExecutionHandle
}
import cromwell.backend.google.pipelines.common.PipelinesApiAsyncBackendJobExecutionActor.JesPendingExecutionHandle
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestManager.PAPIStatusPollRequest
import cromwell.backend.google.pipelines.common.api.RunStatus.UnsuccessfulRunStatus
import cromwell.backend.google.pipelines.common.io.{DiskType, PipelinesApiWorkingDisk}
import cromwell.backend.io.JobPathsSpecHelper._
import cromwell.backend.standard.{
  DefaultStandardAsyncExecutionActorParams,
  StandardAsyncExecutionActorParams,
  StandardAsyncJob,
  StandardExpressionFunctionsParams
}
import cromwell.core._
import cromwell.core.callcaching.NoDocker
import cromwell.core.labels.Labels
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.filesystems.gcs.{GcsPath, GcsPathBuilder, MockGcsPathBuilder}
import cromwell.services.instrumentation.{CromwellBucket, CromwellIncrement}
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.services.keyvalue.InMemoryKvServiceActor
import cromwell.services.keyvalue.KeyValueServiceActor.{KvGet, KvJobKey, KvPair, ScopedKey}
import cromwell.util.JsonFormatting.WomValueJsonFormatter._
import cromwell.util.SampleWdl
import org.scalatest._
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.Tables.Table
import org.slf4j.Logger
import spray.json._
import wdl.transforms.draft2.wdlom2wom.WdlDraft2WomExecutableMakers._
import wdl.transforms.draft2.wdlom2wom._
import wom.WomFileMapper
import wom.expression.NoIoFunctionSet
import wom.graph.CommandCallNode
import wom.transforms.WomExecutableMaker.ops._
import wom.transforms.WomWorkflowDefinitionMaker.ops._
import wom.types._
import wom.values._

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext, Future, Promise}
import scala.language.postfixOps
import scala.util.Success

class PipelinesApiAsyncBackendJobExecutionActorSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with ImplicitSender
    with BackendSpec
    with BeforeAndAfter
    with DefaultJsonProtocol {

  val mockPathBuilder: GcsPathBuilder = MockGcsPathBuilder.instance
  import MockGcsPathBuilder._
  var kvService: ActorRef = system.actorOf(Props(new InMemoryKvServiceActor), "kvService")

  private def gcsPath(str: String) = mockPathBuilder.build(str).getOrElse(fail(s"Invalid gcs path: $str"))

  import PipelinesApiTestConfig._

  implicit val Timeout: FiniteDuration = 25.seconds.dilated

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
       |    docker: "ubuntu:latest"
       |    [PREEMPTIBLE]
       |  }
       |}
       |
       |workflow wf_sup {
       |  call sup
       |}
    """.stripMargin

  val Inputs: Map[FullyQualifiedName, WomValue] = Map("wf_sup.sup.addressee" -> WomString("dog"))

  private val NoOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))

  private lazy val TestableCallContext =
    CallContext(mockPathBuilder.build("gs://root").get, DummyStandardPaths, isDocker = false)

  private lazy val TestableStandardExpressionFunctionsParams: StandardExpressionFunctionsParams =
    new StandardExpressionFunctionsParams {
      override lazy val pathBuilders: List[PathBuilder] = List(mockPathBuilder)
      override lazy val callContext: CallContext = TestableCallContext
      override val ioActorProxy: ActorRef = simpleIoActor
      override val executionContext: ExecutionContext = system.dispatcher
    }

  lazy val TestableJesExpressionFunctions: PipelinesApiExpressionFunctions =
    new PipelinesApiExpressionFunctions(TestableStandardExpressionFunctionsParams)

  private def buildInitializationData(jobDescriptor: BackendJobDescriptor, configuration: PipelinesApiConfiguration) = {
    val pathBuilders =
      Await.result(configuration.configurationDescriptor.pathBuilders(WorkflowOptions.empty), 5.seconds)
    val workflowPaths = PipelinesApiWorkflowPaths(
      jobDescriptor.workflowDescriptor,
      NoCredentials.getInstance(),
      NoCredentials.getInstance(),
      configuration,
      pathBuilders,
      PipelinesApiInitializationActor.defaultStandardStreamNameToFileNameMetadataMapper
    )
    val runtimeAttributesBuilder = PipelinesApiRuntimeAttributes.runtimeAttributesBuilder(configuration)
    val requestFactory: PipelinesApiRequestFactory = new PipelinesApiRequestFactory {
      override def cancelRequest(job: StandardAsyncJob): HttpRequest = null
      override def getRequest(job: StandardAsyncJob): HttpRequest = null
      override def runRequest(createPipelineParameters: PipelinesApiRequestFactory.CreatePipelineParameters,
                              jobLogger: JobLogger
      ): HttpRequest = null
    }
    PipelinesApiBackendInitializationData(workflowPaths,
                                          runtimeAttributesBuilder,
                                          configuration,
                                          null,
                                          requestFactory,
                                          None,
                                          None,
                                          None
    )
  }

  class TestablePipelinesApiJobExecutionActor(params: StandardAsyncExecutionActorParams,
                                              functions: PipelinesApiExpressionFunctions
  ) extends PipelinesApiAsyncBackendJobExecutionActor(params) {

    def this(jobDescriptor: BackendJobDescriptor,
             promise: Promise[BackendJobExecutionResponse],
             jesConfiguration: PipelinesApiConfiguration,
             functions: PipelinesApiExpressionFunctions = TestableJesExpressionFunctions,
             jesSingletonActor: ActorRef = emptyActor,
             ioActor: ActorRef = mockIoActor,
             serviceRegistryActor: ActorRef = kvService
    ) =
      this(
        DefaultStandardAsyncExecutionActorParams(
          jobIdKey = PipelinesApiAsyncBackendJobExecutionActor.JesOperationIdKey,
          serviceRegistryActor = serviceRegistryActor,
          ioActor = ioActor,
          jobDescriptor = jobDescriptor,
          configurationDescriptor = jesConfiguration.configurationDescriptor,
          backendInitializationDataOption = Option(buildInitializationData(jobDescriptor, jesConfiguration)),
          backendSingletonActorOption = Option(jesSingletonActor),
          completionPromise = promise,
          minimumRuntimeSettings = MinimumRuntimeSettings()
        ),
        functions
      )

    override lazy val jobLogger: JobLogger = new JobLogger(
      loggerName = "TestLogger",
      workflowIdForLogging = workflowId.toPossiblyNotRoot,
      rootWorkflowIdForLogging = workflowId.toRoot,
      jobTag = jobTag,
      akkaLogger = Option(log)
    ) {
      override def tag: String = s"$name [UUID(${workflowId.shortString})$jobTag]"
      override val slf4jLoggers: Set[Logger] = Set.empty
    }

    override lazy val backendEngineFunctions: PipelinesApiExpressionFunctions = functions
  }

  private val runtimeAttributesBuilder = PipelinesApiRuntimeAttributes.runtimeAttributesBuilder(papiConfiguration)
  private val workingDisk = PipelinesApiWorkingDisk(DiskType.SSD, 200)

  val DockerAndDiskRuntime: String =
    """
      |runtime {
      |  docker: "ubuntu:latest"
      |  disks: "local-disk 200 SSD"
      |}
    """.stripMargin

  private def buildPreemptibleJobDescriptor(preemptible: Int,
                                            previousPreemptions: Int,
                                            previousUnexpectedRetries: Int,
                                            failedRetriesCountOpt: Option[Int] = None
  ): BackendJobDescriptor = {
    val attempt = previousPreemptions + previousUnexpectedRetries + 1
    val wdlNamespace = WdlNamespaceWithWorkflow
      .load(YoSup.replace("[PREEMPTIBLE]", s"preemptible: $preemptible"), Seq.empty[Draft2ImportResolver])
      .get
    val womDefinition = wdlNamespace.workflow
      .toWomWorkflowDefinition(isASubworkflow = false)
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
        val prefetchedKvEntries = Map(
          PipelinesApiBackendLifecycleActorFactory.preemptionCountKey -> KvPair(
            ScopedKey(workflowDescriptor.id,
                      KvJobKey(key),
                      PipelinesApiBackendLifecycleActorFactory.preemptionCountKey
            ),
            previousPreemptions.toString
          ),
          PipelinesApiBackendLifecycleActorFactory.unexpectedRetryCountKey -> KvPair(
            ScopedKey(workflowDescriptor.id,
                      KvJobKey(key),
                      PipelinesApiBackendLifecycleActorFactory.unexpectedRetryCountKey
            ),
            previousUnexpectedRetries.toString
          )
        )
        val prefetchedKvEntriesUpd = if (failedRetriesCountOpt.isEmpty) {
          prefetchedKvEntries
        } else {
          prefetchedKvEntries + (BackendLifecycleActorFactory.FailedRetryCountKey -> KvPair(
            ScopedKey(workflowDescriptor.id, KvJobKey(key), BackendLifecycleActorFactory.FailedRetryCountKey),
            failedRetriesCountOpt.get.toString
          ))
        }
        BackendJobDescriptor(workflowDescriptor,
                             key,
                             runtimeAttributes,
                             fqnWdlMapToDeclarationMap(Inputs),
                             NoDocker,
                             None,
                             prefetchedKvEntriesUpd
        )
      case Left(badtimes) => fail(badtimes.toList.mkString(", "))
    }
  }

  private case class DockerImageCacheTestingParameters(dockerImageCacheDiskOpt: Option[String],
                                                       dockerImageAsSpecifiedByUser: String,
                                                       isDockerImageCacheUsageRequested: Boolean
  )

  private def executionActor(jobDescriptor: BackendJobDescriptor,
                             promise: Promise[BackendJobExecutionResponse],
                             jesSingletonActor: ActorRef,
                             shouldBePreemptible: Boolean,
                             serviceRegistryActor: ActorRef = kvService,
                             referenceInputFilesOpt: Option[Set[PipelinesApiInput]] = None,
                             dockerImageCacheTestingParamsOpt: Option[DockerImageCacheTestingParameters] = None
  ): ActorRef = {

    val job = StandardAsyncJob(UUID.randomUUID().toString)
    val run = Run(job)
    val handle = new JesPendingExecutionHandle(jobDescriptor, run.job, Option(run), None)

    class ExecuteOrRecoverActor
        extends TestablePipelinesApiJobExecutionActor(jobDescriptor,
                                                      promise,
                                                      papiConfiguration,
                                                      jesSingletonActor = jesSingletonActor,
                                                      serviceRegistryActor = serviceRegistryActor
        ) {
      override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
        sendIncrementMetricsForReferenceFiles(referenceInputFilesOpt)
        dockerImageCacheTestingParamsOpt.foreach { dockerImageCacheTestingParams =>
          sendIncrementMetricsForDockerImageCache(
            dockerImageCacheTestingParams.dockerImageCacheDiskOpt,
            dockerImageCacheTestingParams.dockerImageAsSpecifiedByUser,
            dockerImageCacheTestingParams.isDockerImageCacheUsageRequested
          )
        }

        if (preemptible == shouldBePreemptible) Future.successful(handle)
        else Future.failed(new Exception(s"Test expected preemptible to be $shouldBePreemptible but got $preemptible"))
      }
    }

    system.actorOf(Props(new ExecuteOrRecoverActor), "ExecuteOrRecoverActor-" + UUID.randomUUID)
  }

  private def runAndFail(previousPreemptions: Int,
                         previousUnexpectedRetries: Int,
                         preemptible: Int,
                         errorCode: Status,
                         innerErrorMessage: String,
                         expectPreemptible: Boolean
  ): BackendJobExecutionResponse = {

    val runStatus = UnsuccessfulRunStatus(errorCode,
                                          Option(innerErrorMessage),
                                          Seq.empty,
                                          Option("fakeMachine"),
                                          Option("fakeZone"),
                                          Option("fakeInstance"),
                                          expectPreemptible
    )
    val statusPoller = TestProbe("statusPoller")

    val promise = Promise[BackendJobExecutionResponse]()
    val jobDescriptor = buildPreemptibleJobDescriptor(preemptible, previousPreemptions, previousUnexpectedRetries)

    // TODO: Use this to check the new KV entries are there!
    // val kvProbe = TestProbe("kvProbe")

    val backend = executionActor(jobDescriptor, promise, statusPoller.ref, expectPreemptible)
    backend ! Execute
    statusPoller.expectMsgPF(max = Timeout, hint = "awaiting status poll") { case _: PAPIStatusPollRequest =>
      backend ! runStatus
    }

    Await.result(promise.future, Timeout)
  }

  def buildPreemptibleTestActorRef(attempt: Int,
                                   preemptible: Int,
                                   failedRetriesCountOpt: Option[Int] = None
  ): TestActorRef[TestablePipelinesApiJobExecutionActor] = {
    // For this test we say that all previous attempts were preempted:
    val jobDescriptor = buildPreemptibleJobDescriptor(preemptible,
                                                      attempt - 1,
                                                      previousUnexpectedRetries = 0,
                                                      failedRetriesCountOpt = failedRetriesCountOpt
    )
    val props = Props(
      new TestablePipelinesApiJobExecutionActor(jobDescriptor,
                                                Promise(),
                                                papiConfiguration,
                                                TestableJesExpressionFunctions,
                                                emptyActor,
                                                failIoActor
      )
    )
    TestActorRef(props, s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")
  }

  behavior of "JesAsyncBackendJobExecutionActor"

  private val timeout = 25 seconds

  { // Set of "handle call failures appropriately with respect to preemption and failure" tests
    val expectations = Table(
      ("previous_preemptions",
       "previous_unexpectedRetries",
       "preemptible",
       "errorCode",
       "message",
       "shouldRunAsPreemptible",
       "shouldRetry"
      ),
      // No preemptible attempts allowed, but standard failures should be retried.
      (0,
       0,
       0,
       Status.ABORTED,
       "13: retryable error",
       false,
       true
      ), // This is the new "unexpected failure" mode, which is now retried
      (0, 1, 0, Status.ABORTED, "13: retryable error", false, true),
      (0, 2, 0, Status.ABORTED, "13: retryable error", false, false), // The third unexpected failure is a real failure.
      (0,
       0,
       0,
       Status.ABORTED,
       "14: usually means preempted...?",
       false,
       false
      ), // Usually means "preempted', but this wasn't a preemptible VM, so this should just be a failure.
      (0, 0, 0, Status.ABORTED, "15: other error", false, false),
      (0, 0, 0, Status.OUT_OF_RANGE, "13: unexpected error", false, false),
      (0, 0, 0, Status.OUT_OF_RANGE, "14: test error msg", false, false),
      // These commented out tests should be uncommented if/when we stop mapping 13 to 14 in preemption mode
      // 1 preemptible attempt allowed, but not all failures represent preemptions.
//      (0, 0, 1, Status.ABORTED, "13: retryable error", true, true),
//      (0, 1, 1, Status.ABORTED, "13: retryable error", true, true),
//      (0, 2, 1, Status.ABORTED, "13: retryable error", true, false),
      // The following 13 based test should be removed if/when we stop mapping 13 to 14 in preemption mode
      (0, 0, 1, Status.ABORTED, "13: retryable error", true, true),
      (0, 0, 1, Status.ABORTED, "14: preempted", true, true),
      (0, 0, 1, Status.UNKNOWN, "Instance failed to start due to preemption.", true, true),
      (0, 0, 1, Status.ABORTED, "15: other error", true, false),
      (0, 0, 1, Status.OUT_OF_RANGE, "13: retryable error", true, false),
      (0, 0, 1, Status.OUT_OF_RANGE, "14: preempted", true, false),
      (0, 0, 1, Status.OUT_OF_RANGE, "Instance failed to start due to preemption.", true, false),
      // 1 preemptible attempt allowed, but since we're now on the second preemption attempt only 13s should be retryable.
      (1, 0, 1, Status.ABORTED, "13: retryable error", false, true),
      (1, 1, 1, Status.ABORTED, "13: retryable error", false, true),
      (1, 2, 1, Status.ABORTED, "13: retryable error", false, false),
      (1, 0, 1, Status.ABORTED, "14: preempted", false, false),
      (1, 0, 1, Status.UNKNOWN, "Instance failed to start due to preemption.", false, false),
      (1, 0, 1, Status.ABORTED, "15: other error", false, false),
      (1, 0, 1, Status.OUT_OF_RANGE, "13: retryable error", false, false),
      (1, 0, 1, Status.OUT_OF_RANGE, "14: preempted", false, false),
      (1, 0, 1, Status.OUT_OF_RANGE, "Instance failed to start due to preemption.", false, false)
    )

    expectations foreach {
      case (previousPreemptions,
            previousUnexpectedRetries,
            preemptible,
            errorCode,
            innerErrorMessage,
            shouldBePreemptible,
            shouldRetry
          ) =>
        val descriptor =
          s"previousPreemptions=$previousPreemptions, previousUnexpectedRetries=$previousUnexpectedRetries preemptible=$preemptible, errorCode=$errorCode, innerErrorMessage=$innerErrorMessage"
        it should s"handle call failures appropriately with respect to preemption and failure ($descriptor)" in {
          runAndFail(previousPreemptions,
                     previousUnexpectedRetries,
                     preemptible,
                     errorCode,
                     innerErrorMessage,
                     shouldBePreemptible
          ) match {
            case response: JobFailedNonRetryableResponse =>
              if (shouldRetry)
                fail(s"A should-be-retried job ($descriptor) was sent back to the engine with: $response")
            case response: JobFailedRetryableResponse =>
              if (!shouldRetry)
                fail(s"A shouldn't-be-retried job ($descriptor) was sent back to the engine with $response")
            case huh => fail(s"Unexpected response: $huh")
          }
        }
    }
  }

  it should "send proper value for \"number of reference files used gauge\" metric, or don't send anything if reference disks feature is disabled" in {
    val expectedInput1 = PipelinesApiFileInput(name = "testfile1",
                                               relativeHostPath =
                                                 DefaultPathBuilder.build(Paths.get(s"test/reference/path/file1")),
                                               mount = null,
                                               cloudPath = null
    )
    val expectedInput2 = PipelinesApiFileInput(name = "testfile2",
                                               relativeHostPath =
                                                 DefaultPathBuilder.build(Paths.get(s"test/reference/path/file2")),
                                               mount = null,
                                               cloudPath = null
    )
    val expectedReferenceInputFiles = Set[PipelinesApiInput](expectedInput1, expectedInput2)

    val expectedMsg1 = InstrumentationServiceMessage(
      CromwellIncrement(
        CromwellBucket(List.empty, NonEmptyList.of("referencefiles", expectedInput1.relativeHostPath.pathAsString))
      )
    )
    val expectedMsg2 = InstrumentationServiceMessage(
      CromwellIncrement(
        CromwellBucket(List.empty, NonEmptyList.of("referencefiles", expectedInput2.relativeHostPath.pathAsString))
      )
    )

    val jobDescriptor = buildPreemptibleJobDescriptor(0, 0, 0)
    val serviceRegistryProbe = TestProbe()

    val backend1 = executionActor(
      jobDescriptor,
      Promise[BackendJobExecutionResponse](),
      TestProbe().ref,
      shouldBePreemptible = false,
      serviceRegistryActor = serviceRegistryProbe.ref,
      referenceInputFilesOpt = Option(expectedReferenceInputFiles)
    )
    backend1 ! Execute
    serviceRegistryProbe.expectMsgAllOf(expectedMsg1, expectedMsg2)

    val backend2 = executionActor(
      jobDescriptor,
      Promise[BackendJobExecutionResponse](),
      TestProbe().ref,
      shouldBePreemptible = false,
      serviceRegistryActor = serviceRegistryProbe.ref,
      referenceInputFilesOpt = None
    )
    backend2 ! Execute
    serviceRegistryProbe.expectNoMessage(timeout)
  }

  it should "sends proper metrics for docker image cache feature" in {
    val jobDescriptor = buildPreemptibleJobDescriptor(0, 0, 0)
    val serviceRegistryProbe = TestProbe()
    val madeUpDockerImageName = "test_madeup_docker_image_name"

    val expectedMessageWhenRequestedNotFound = InstrumentationServiceMessage(
      CromwellIncrement(
        CromwellBucket(List.empty,
                       NonEmptyList("docker", List("image", "cache", "image_not_in_cache", madeUpDockerImageName))
        )
      )
    )
    val backendDockerCacheRequestedButNotFound = executionActor(
      jobDescriptor,
      Promise[BackendJobExecutionResponse](),
      TestProbe().ref,
      shouldBePreemptible = false,
      serviceRegistryActor = serviceRegistryProbe.ref,
      dockerImageCacheTestingParamsOpt = Option(
        DockerImageCacheTestingParameters(
          None,
          "test_madeup_docker_image_name",
          isDockerImageCacheUsageRequested = true
        )
      )
    )
    backendDockerCacheRequestedButNotFound ! Execute
    serviceRegistryProbe.expectMsg(expectedMessageWhenRequestedNotFound)

    val expectedMessageWhenRequestedAndFound = InstrumentationServiceMessage(
      CromwellIncrement(
        CromwellBucket(List.empty,
                       NonEmptyList("docker", List("image", "cache", "used_image_from_cache", madeUpDockerImageName))
        )
      )
    )
    val backendDockerCacheRequestedAndFound = executionActor(
      jobDescriptor,
      Promise[BackendJobExecutionResponse](),
      TestProbe().ref,
      shouldBePreemptible = false,
      serviceRegistryActor = serviceRegistryProbe.ref,
      dockerImageCacheTestingParamsOpt = Option(
        DockerImageCacheTestingParameters(
          Option("test_madeup_disk_image_name"),
          "test_madeup_docker_image_name",
          isDockerImageCacheUsageRequested = true
        )
      )
    )
    backendDockerCacheRequestedAndFound ! Execute
    serviceRegistryProbe.expectMsg(expectedMessageWhenRequestedAndFound)

    val expectedMessageWhenNotRequestedButFound = InstrumentationServiceMessage(
      CromwellIncrement(
        CromwellBucket(List.empty,
                       NonEmptyList("docker", List("image", "cache", "cached_image_not_used", madeUpDockerImageName))
        )
      )
    )
    val backendDockerCacheNotRequestedButFound = executionActor(
      jobDescriptor,
      Promise[BackendJobExecutionResponse](),
      TestProbe().ref,
      shouldBePreemptible = false,
      serviceRegistryActor = serviceRegistryProbe.ref,
      dockerImageCacheTestingParamsOpt = Option(
        DockerImageCacheTestingParameters(
          Option("test_madeup_disk_image_name"),
          "test_madeup_docker_image_name",
          isDockerImageCacheUsageRequested = false
        )
      )
    )
    backendDockerCacheNotRequestedButFound ! Execute
    serviceRegistryProbe.expectMsg(expectedMessageWhenNotRequestedButFound)

    val backendDockerCacheNotRequestedNotFound = executionActor(
      jobDescriptor,
      Promise[BackendJobExecutionResponse](),
      TestProbe().ref,
      shouldBePreemptible = false,
      serviceRegistryActor = serviceRegistryProbe.ref,
      dockerImageCacheTestingParamsOpt = Option(
        DockerImageCacheTestingParameters(
          None,
          "test_madeup_docker_image_name",
          isDockerImageCacheUsageRequested = false
        )
      )
    )
    backendDockerCacheNotRequestedNotFound ! Execute
    serviceRegistryProbe.expectNoMessage(timeout)
  }

  it should "not restart 2 of 1 unexpected shutdowns without another preemptible VM" in {
    val actorRef = buildPreemptibleTestActorRef(2, 1)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(
      Status.ABORTED,
      Option("14: VM XXX shut down unexpectedly."),
      Seq.empty,
      Option("fakeMachine"),
      Option("fakeZone"),
      Option("fakeInstance"),
      wasPreemptible = true
    )
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    val failedHandle = result.asInstanceOf[FailedNonRetryableExecutionHandle]
    failedHandle.returnCode shouldBe None
  }

  it should "restart 1 of 1 unexpected shutdowns without another preemptible VM" in {
    val actorRef = buildPreemptibleTestActorRef(1, 1)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(
      Status.ABORTED,
      Option("14: VM XXX shut down unexpectedly."),
      Seq.empty,
      Option("fakeMachine"),
      Option("fakeZone"),
      Option("fakeInstance"),
      wasPreemptible = true
    )
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    retryableHandle.returnCode shouldBe None
    retryableHandle.throwable.getMessage should include("will be restarted with a non-preemptible VM")
  }

  it should "restart 1 of 2 unexpected shutdowns with another preemptible VM" in {
    val actorRef = buildPreemptibleTestActorRef(1, 2)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(
      Status.ABORTED,
      Option("14: VM XXX shut down unexpectedly."),
      Seq.empty,
      Option("fakeMachine"),
      Option("fakeZone"),
      Option("fakeInstance"),
      wasPreemptible = true
    )
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    retryableHandle.returnCode shouldBe None
    retryableHandle.throwable.getMessage should include("will be restarted with another preemptible VM")
  }

  it should "treat a JES message 13 as preemptible if the VM was preemptible" in {
    val actorRef = buildPreemptibleTestActorRef(1, 2)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(Status.ABORTED,
                                             Option("13: Retryable Error."),
                                             Seq.empty,
                                             Option("fakeMachine"),
                                             Option("fakeZone"),
                                             Option("fakeInstance"),
                                             wasPreemptible = true
    )
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    retryableHandle.returnCode shouldBe None
    retryableHandle.throwable.getMessage should include("will be restarted with another preemptible VM")
  }

  it should "treat a PAPI v2 style error message as preemptible if the VM was preemptible" in {
    val actorRef = buildPreemptibleTestActorRef(1, 2)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(
      Status.ABORTED,
      Option(PipelinesApiAsyncBackendJobExecutionActor.FailedV2Style),
      Seq.empty,
      Option("fakeMachine"),
      Option("fakeZone"),
      Option("fakeInstance"),
      wasPreemptible = true
    )
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    retryableHandle.returnCode shouldBe None
    retryableHandle.throwable.getMessage should include("will be restarted with another preemptible VM")
  }

  it should "when at the preemptible limit restart a PAPI v2 style error message with a non-preemptible VM" in {
    val actorRef = buildPreemptibleTestActorRef(1, 1)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(
      Status.ABORTED,
      Option(PipelinesApiAsyncBackendJobExecutionActor.FailedV2Style),
      Seq.empty,
      Option("fakeMachine"),
      Option("fakeZone"),
      Option("fakeInstance"),
      wasPreemptible = true
    )
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
    val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
    retryableHandle.returnCode shouldBe None
    retryableHandle.throwable.getMessage should include("The call will be restarted with a non-preemptible VM.")
  }

  it should "merge KvPairs from previous attempt with KvPairs generated by backend for the next attempt and send to KeyValueServiceActor in case of job retryable failure" in {
    val actorRef = buildPreemptibleTestActorRef(attempt = 1, preemptible = 1, failedRetriesCountOpt = Option(1))
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    val failedStatus = UnsuccessfulRunStatus(
      Status.ABORTED,
      Option(PipelinesApiAsyncBackendJobExecutionActor.FailedV2Style),
      Seq.empty,
      Option("fakeMachine"),
      Option("fakeZone"),
      Option("fakeInstance"),
      wasPreemptible = true
    )
    val executionResult = jesBackend.handleExecutionResult(failedStatus, handle)
    val result = Await.result(executionResult, timeout)
    result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true

    val probe = TestProbe("probe")
    val job = jesBackend.workflowDescriptor.callable.taskCallNodes.head
    val key = BackendJobDescriptorKey(job, None, attempt = 2)

    val scopedKeyPreempt =
      ScopedKey(jesBackend.workflowId, KvJobKey(key), PipelinesApiBackendLifecycleActorFactory.preemptionCountKey)
    probe.send(kvService, KvGet(scopedKeyPreempt))
    probe.expectMsgPF(5.seconds) { case kvPreempt: KvPair =>
      kvPreempt.value shouldBe 1.toString
    }

    val scopedKeyUnexpected =
      ScopedKey(jesBackend.workflowId, KvJobKey(key), PipelinesApiBackendLifecycleActorFactory.unexpectedRetryCountKey)
    probe.send(kvService, KvGet(scopedKeyUnexpected))
    probe.expectMsgPF(5.seconds) { case kvUnexpected: KvPair =>
      kvUnexpected.value shouldBe 0.toString
    }

    val scopedKeyFailedAttempts =
      ScopedKey(jesBackend.workflowId, KvJobKey(key), BackendLifecycleActorFactory.FailedRetryCountKey)
    probe.send(kvService, KvGet(scopedKeyFailedAttempts))
    probe.expectMsgPF(5.seconds) { case kvFailedRetry: KvPair =>
      kvFailedRetry.value shouldBe 1.toString
    }
  }

  it should "handle Failure Status for various errors" in {
    val actorRef = buildPreemptibleTestActorRef(1, 1)
    val jesBackend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new JesPendingExecutionHandle(null, runId, None, None)

    def checkFailedResult(errorCode: Status, errorMessage: Option[String]): ExecutionHandle = {
      val failed = UnsuccessfulRunStatus(errorCode,
                                         errorMessage,
                                         Seq.empty,
                                         Option("fakeMachine"),
                                         Option("fakeZone"),
                                         Option("fakeInstance"),
                                         wasPreemptible = true
      )
      Await.result(jesBackend.handleExecutionResult(failed, handle), timeout)
    }

    checkFailedResult(Status.ABORTED, Option("15: Other type of error."))
      .isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Status.OUT_OF_RANGE, Option("14: Wrong errorCode."))
      .isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Status.ABORTED, Option("Weird error message."))
      .isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Status.ABORTED, Option("UnparsableInt: Even weirder error message."))
      .isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Status.ABORTED, None).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Status.CANCELLED, Option("Operation canceled at")) shouldBe AbortedExecutionHandle

    actorRef.stop()
  }

  it should "map GCS paths and *only* GCS paths to local" in {
    val wdlString =
      s"""|workflow wf {
          |  call t
          |}
          |
          |task t {
          |  String abc
          |  File lf
          |  File gcsf
          |  command {}
          |  runtime { docker: "ubuntu" }
          |}
          |""".stripMargin

    val stringKey = "wf.t.abc"
    val stringVal = WomString("abc")
    val localFileKey = "wf.t.lf"
    val localFileVal = WomSingleFile("/blah/abc")
    val gcsFileKey = "wf.t.gcsf"
    val gcsFileVal = WomSingleFile("gs://blah/abc")

    val inputs: Map[String, WomValue] = Map(
      stringKey -> stringVal,
      localFileKey -> localFileVal,
      gcsFileKey -> gcsFileVal
    )

    val wdlNamespace = WdlNamespaceWithWorkflow
      .load(
        wdlString,
        Seq.empty[Draft2ImportResolver]
      )
      .get
    val womWorkflow = wdlNamespace.workflow
      .toWomWorkflowDefinition(isASubworkflow = false)
      .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    wdlNamespace.toWomExecutable(Option(inputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap { case (port, v) =>
          v.select[WomValue] map { port -> _ }
        }

        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val call: CommandCallNode = workflowDescriptor.callable.graph.nodes.collectFirst { case t: CommandCallNode =>
          t
        }.get
        val key = BackendJobDescriptorKey(call, None, 1)
        val runtimeAttributes = makeRuntimeAttributes(call)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor,
                                                 key,
                                                 runtimeAttributes,
                                                 fqnWdlMapToDeclarationMap(inputs),
                                                 NoDocker,
                                                 None,
                                                 Map.empty
        )

        val props = Props(new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration))
        val testActorRef = TestActorRef[TestablePipelinesApiJobExecutionActor](
          props,
          s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
        )

        def gcsPathToLocal(womValue: WomValue): WomValue =
          WomFileMapper.mapWomFiles(testActorRef.underlyingActor.mapCommandLineWomFile)(womValue).get

        val mappedInputs = jobDescriptor.localInputs safeMapValues gcsPathToLocal

        mappedInputs(stringKey) match {
          case WomString(v) => assert(v.equalsIgnoreCase(stringVal.value))
          case _ => fail("test setup error")
        }

        mappedInputs(localFileKey) match {
          case wdlFile: WomSingleFile => assert(wdlFile.value.equalsIgnoreCase(localFileVal.value))
          case _ => fail("test setup error")
        }

        mappedInputs(gcsFileKey) match {
          case wdlFile: WomSingleFile => assert(wdlFile.value.equalsIgnoreCase("/cromwell_root/blah/abc"))
          case _ => fail("test setup error")
        }
      case Left(badtimes) => fail(badtimes.toList.mkString(", "))
    }
  }

  private val dockerAndDiskMapsWdlNamespace =
    WdlNamespaceWithWorkflow
      .load(
        SampleWdl.CurrentDirectoryMaps.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
        Seq.empty[Draft2ImportResolver]
      )
      .get

  private val dockerAndDiskArrayWdlNamespace =
    WdlNamespaceWithWorkflow
      .load(
        SampleWdl.CurrentDirectoryArray.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
        Seq.empty[Draft2ImportResolver]
      )
      .get

  private val dockerAndDiskFilesWdlNamespace =
    WdlNamespaceWithWorkflow
      .load(
        SampleWdl.CurrentDirectoryFiles.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
        Seq.empty[Draft2ImportResolver]
      )
      .get

  it should "generate correct JesFileInputs from a WdlMap" in {
    val inputs: Map[String, WomValue] = Map(
      "stringToFileMap" -> WomMap(
        WomMapType(WomStringType, WomSingleFileType),
        Map(
          WomString("stringTofile1") -> WomSingleFile("gs://path/to/stringTofile1"),
          WomString("stringTofile2") -> WomSingleFile("gs://path/to/stringTofile2")
        )
      ),
      "fileToStringMap" -> WomMap(
        WomMapType(WomSingleFileType, WomStringType),
        Map(
          WomSingleFile("gs://path/to/fileToString1") -> WomString("fileToString1"),
          WomSingleFile("gs://path/to/fileToString2") -> WomString("fileToString2")
        )
      ),
      "fileToFileMap" -> WomMap(
        WomMapType(WomSingleFileType, WomSingleFileType),
        Map(
          WomSingleFile("gs://path/to/fileToFile1Key") -> WomSingleFile("gs://path/to/fileToFile1Value"),
          WomSingleFile("gs://path/to/fileToFile2Key") -> WomSingleFile("gs://path/to/fileToFile2Value")
        )
      ),
      "stringToString" -> WomMap(
        WomMapType(WomStringType, WomStringType),
        Map(
          WomString("stringToString1") -> WomString("path/to/stringToString1"),
          WomString("stringToString2") -> WomString("path/to/stringToString2")
        )
      )
    )

    val workflowInputs = inputs map { case (key, value) =>
      (s"wf_whereami.whereami.$key", value)
    }

    val womWorkflow =
      dockerAndDiskMapsWdlNamespace.workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    val womExecutableChecked =
      dockerAndDiskMapsWdlNamespace
        .toWomExecutable(Option(workflowInputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true)
    womExecutableChecked match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap { case (port, v) =>
          v.select[WomValue] map { port -> _ }
        }
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor,
                                                 key,
                                                 runtimeAttributes,
                                                 fqnWdlMapToDeclarationMap(inputs),
                                                 NoDocker,
                                                 None,
                                                 Map.empty
        )

        val props = Props(new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration))
        val testActorRef = TestActorRef[TestablePipelinesApiJobExecutionActor](
          props,
          s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
        )

        val jesInputs = testActorRef.underlyingActor.generateInputs(jobDescriptor)
        jesInputs should have size 8
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.stringToFileMap-0",
            cloudPath = gcsPath("gs://path/to/stringTofile1"),
            relativeHostPath = DefaultPathBuilder.get("path/to/stringTofile1"),
            mount = workingDisk
          )
        )
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.stringToFileMap-1",
            cloudPath = gcsPath("gs://path/to/stringTofile2"),
            relativeHostPath = DefaultPathBuilder.get("path/to/stringTofile2"),
            mount = workingDisk
          )
        )
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.fileToStringMap-0",
            cloudPath = gcsPath("gs://path/to/fileToString1"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToString1"),
            mount = workingDisk
          )
        )
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.fileToStringMap-1",
            cloudPath = gcsPath("gs://path/to/fileToString2"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToString2"),
            mount = workingDisk
          )
        )
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.fileToFileMap-0",
            cloudPath = gcsPath("gs://path/to/fileToFile1Key"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToFile1Key"),
            mount = workingDisk
          )
        )
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.fileToFileMap-1",
            cloudPath = gcsPath("gs://path/to/fileToFile1Value"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToFile1Value"),
            mount = workingDisk
          )
        )
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.fileToFileMap-2",
            cloudPath = gcsPath("gs://path/to/fileToFile2Key"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToFile2Key"),
            mount = workingDisk
          )
        )
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.fileToFileMap-3",
            cloudPath = gcsPath("gs://path/to/fileToFile2Value"),
            relativeHostPath = DefaultPathBuilder.get("path/to/fileToFile2Value"),
            mount = workingDisk
          )
        )

      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  private def makeJesActorRef(sampleWdl: SampleWdl,
                              workflowInputs: Map[FullyQualifiedName, WomValue],
                              callName: LocallyQualifiedName,
                              callInputs: Map[LocallyQualifiedName, WomValue],
                              functions: PipelinesApiExpressionFunctions = TestableJesExpressionFunctions
  ): TestActorRef[TestablePipelinesApiJobExecutionActor] = {
    val wdlNamespaceWithWorkflow =
      WdlNamespaceWithWorkflow
        .load(
          sampleWdl.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
          Seq.empty[Draft2ImportResolver]
        )
        .get
    val womWorkflow =
      wdlNamespaceWithWorkflow.workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    val womExecutableChecked =
      wdlNamespaceWithWorkflow
        .toWomExecutable(
          Option(workflowInputs.toJson.compactPrint),
          NoIoFunctionSet,
          strictValidation = true
        )
    womExecutableChecked match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap { case (port, v) =>
          v.select[WomValue] map { port -> _ }
        }
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == callName).get
        val key = BackendJobDescriptorKey(call, None, 1)
        val runtimeAttributes = makeRuntimeAttributes(call)
        val jobDescriptor =
          BackendJobDescriptor(
            workflowDescriptor = workflowDescriptor,
            key = key,
            runtimeAttributes = runtimeAttributes,
            evaluatedTaskInputs = fqnWdlMapToDeclarationMap(callInputs),
            maybeCallCachingEligible = NoDocker,
            dockerSize = None,
            prefetchedKvStoreEntries = Map.empty
          )

        val props = Props(
          new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration, functions)
        )
        TestActorRef[TestablePipelinesApiJobExecutionActor](
          props,
          s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
        )
      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "generate correct JesOutputs" in {
    val womFile = WomSingleFile("gs://blah/b/c.txt")
    val workflowInputs = Map("file_passing.f" -> womFile)
    val callInputs = Map(
      "in" -> womFile, // how does one programmatically map the wf inputs to the call inputs?
      "out_name" -> WomString("out") // is it expected that this isn't using the default?
    )
    val jesBackend = makeJesActorRef(SampleWdl.FilePassingWorkflow, workflowInputs, "a", callInputs).underlyingActor
    val jobDescriptor = jesBackend.jobDescriptor
    val workflowId = jesBackend.workflowId
    val jesInputs = jesBackend.generateInputs(jobDescriptor)
    jesInputs should have size 1
    jesInputs should contain(
      PipelinesApiFileInput(
        name = "file_passing.a.in-0",
        cloudPath = gcsPath("gs://blah/b/c.txt"),
        relativeHostPath = DefaultPathBuilder.get("blah/b/c.txt"),
        mount = workingDisk
      )
    )
    val jesOutputs = jesBackend.generateOutputs(jobDescriptor)
    jesOutputs should have size 1
    jesOutputs should contain(
      PipelinesApiFileOutput(
        "out",
        gcsPath(s"gs://my-cromwell-workflows-bucket/file_passing/$workflowId/call-a/out"),
        DefaultPathBuilder.get("out"),
        workingDisk,
        optional = false,
        secondary = false
      )
    )
  }

  it should "generate correct JesInputs when a command line contains a write_lines call in it" in {
    val inputs = Map(
      "strs" -> WomArray(WomArrayType(WomStringType), Seq("A", "B", "C").map(WomString))
    )

    class TestPipelinesApiExpressionFunctions
        extends PipelinesApiExpressionFunctions(TestableStandardExpressionFunctionsParams) {
      override def writeFile(path: String, content: String): Future[WomSingleFile] =
        Future.fromTry(Success(WomSingleFile(s"gs://some/path/file.txt")))
    }

    val functions = new TestPipelinesApiExpressionFunctions
    val jesBackend = makeJesActorRef(SampleWdl.ArrayIO, Map.empty, "serialize", inputs, functions).underlyingActor
    val jobDescriptor = jesBackend.jobDescriptor
    val jesInputs = jesBackend.generateInputs(jobDescriptor)
    jesInputs should have size 1
    jesInputs should contain(
      PipelinesApiFileInput(
        name = "c35ad8d3-0",
        cloudPath = gcsPath("gs://some/path/file.txt"),
        relativeHostPath = DefaultPathBuilder.get("some/path/file.txt"),
        mount = workingDisk
      )
    )
    val jesOutputs = jesBackend.generateOutputs(jobDescriptor)
    jesOutputs should have size 0
  }

  it should "generate correct JesFileInputs from a WdlArray" in {
    val inputs: Map[String, WomValue] = Map(
      "fileArray" ->
        WomArray(WomArrayType(WomSingleFileType),
                 Seq(WomSingleFile("gs://path/to/file1"), WomSingleFile("gs://path/to/file2"))
        )
    )

    val workflowInputs = inputs map { case (key, value) =>
      (s"wf_whereami.whereami.$key", value)
    }

    val womWorkflow =
      dockerAndDiskArrayWdlNamespace.workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    val womExecutableChecked =
      dockerAndDiskArrayWdlNamespace
        .toWomExecutable(Option(workflowInputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true)
    womExecutableChecked match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap { case (port, v) =>
          v.select[WomValue] map { port -> _ }
        }
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor,
                                                 key,
                                                 runtimeAttributes,
                                                 fqnWdlMapToDeclarationMap(inputs),
                                                 NoDocker,
                                                 None,
                                                 Map.empty
        )

        val props = Props(new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration))
        val testActorRef = TestActorRef[TestablePipelinesApiJobExecutionActor](
          props,
          s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
        )

        val jesInputs = testActorRef.underlyingActor.generateInputs(jobDescriptor)
        jesInputs should have size 2
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.fileArray-0",
            cloudPath = gcsPath("gs://path/to/file1"),
            relativeHostPath = DefaultPathBuilder.get("path/to/file1"),
            mount = workingDisk
          )
        )
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.fileArray-1",
            cloudPath = gcsPath("gs://path/to/file2"),
            relativeHostPath = DefaultPathBuilder.get("path/to/file2"),
            mount = workingDisk
          )
        )
      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "generate correct JesFileInputs from a WdlFile" in {
    val inputs: Map[String, WomValue] = Map(
      "file1" -> WomSingleFile("gs://path/to/file1"),
      "file2" -> WomSingleFile("gs://path/to/file2")
    )

    val workflowInputs = inputs map { case (key, value) =>
      (s"wf_whereami.whereami.$key", value)
    }

    val womWorkflow =
      dockerAndDiskFilesWdlNamespace.workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    val womExecutableChecked =
      dockerAndDiskFilesWdlNamespace
        .toWomExecutable(Option(workflowInputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true)
    womExecutableChecked match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap { case (port, v) =>
          v.select[WomValue] map { port -> _ }
        }
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty,
          HogGroup("foo"),
          List.empty,
          None
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor,
                                                 key,
                                                 runtimeAttributes,
                                                 fqnWdlMapToDeclarationMap(inputs),
                                                 NoDocker,
                                                 None,
                                                 Map.empty
        )

        val props = Props(new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration))
        val testActorRef = TestActorRef[TestablePipelinesApiJobExecutionActor](
          props,
          s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
        )

        val jesInputs = testActorRef.underlyingActor.generateInputs(jobDescriptor)
        jesInputs should have size 2
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.file1-0",
            cloudPath = gcsPath("gs://path/to/file1"),
            relativeHostPath = DefaultPathBuilder.get("path/to/file1"),
            mount = workingDisk
          )
        )
        jesInputs should contain(
          PipelinesApiFileInput(
            name = "wf_whereami.whereami.file2-0",
            cloudPath = gcsPath("gs://path/to/file2"),
            relativeHostPath = DefaultPathBuilder.get("path/to/file2"),
            mount = workingDisk
          )
        )

      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "convert local Paths back to corresponding GCS paths in JesOutputs" in {
    val jesOutputs = Set(
      PipelinesApiFileOutput(
        "/cromwell_root/path/to/file1",
        gcsPath("gs://path/to/file1"),
        DefaultPathBuilder.get("/cromwell_root/path/to/file1"),
        workingDisk,
        optional = false,
        secondary = false
      ),
      PipelinesApiFileOutput(
        "/cromwell_root/path/to/file2",
        gcsPath("gs://path/to/file2"),
        DefaultPathBuilder.get("/cromwell_root/path/to/file2"),
        workingDisk,
        optional = false,
        secondary = false
      ),
      PipelinesApiFileOutput(
        "/cromwell_root/path/to/file3",
        gcsPath("gs://path/to/file3"),
        DefaultPathBuilder.get("/cromwell_root/path/to/file3"),
        workingDisk,
        optional = false,
        secondary = false
      ),
      PipelinesApiFileOutput(
        "/cromwell_root/path/to/file4",
        gcsPath("gs://path/to/file4"),
        DefaultPathBuilder.get("/cromwell_root/path/to/file4"),
        workingDisk,
        optional = false,
        secondary = false
      ),
      PipelinesApiFileOutput(
        "/cromwell_root/path/to/file5",
        gcsPath("gs://path/to/file5"),
        DefaultPathBuilder.get("/cromwell_root/path/to/file5"),
        workingDisk,
        optional = false,
        secondary = false
      )
    )
    val outputValues = Seq(
      WomSingleFile("/cromwell_root/path/to/file1"),
      WomArray(WomArrayType(WomSingleFileType),
               Seq(WomSingleFile("/cromwell_root/path/to/file2"), WomSingleFile("/cromwell_root/path/to/file3"))
      ),
      WomMap(WomMapType(WomSingleFileType, WomSingleFileType),
             Map(
               WomSingleFile("/cromwell_root/path/to/file4") -> WomSingleFile("/cromwell_root/path/to/file5")
             )
      )
    )

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow
        .load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
              Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      NoOptions,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration))
    val testActorRef = TestActorRef[TestablePipelinesApiJobExecutionActor](
      props,
      s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    def wdlValueToGcsPath(jesOutputs: Set[PipelinesApiFileOutput])(womValue: WomValue): WomValue =
      WomFileMapper.mapWomFiles(testActorRef.underlyingActor.womFileToGcsPath(jesOutputs.toSet))(womValue).get

    val result = outputValues map wdlValueToGcsPath(jesOutputs)
    result should have size 3
    result should contain(WomSingleFile("gs://path/to/file1"))
    result should contain(
      WomArray(WomArrayType(WomSingleFileType),
               Seq(WomSingleFile("gs://path/to/file2"), WomSingleFile("gs://path/to/file3"))
      )
    )
    result should contain(
      WomMap(WomMapType(WomSingleFileType, WomSingleFileType),
             Map(WomSingleFile("gs://path/to/file4") -> WomSingleFile("gs://path/to/file5"))
      )
    )
  }

  it should "create a JesFileInput for the monitoring script, when specified" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow
        .load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
              Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions.fromJsonString("""{"monitoring_script": "gs://path/to/script"}""").get,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
    val runtimeAttributes = makeRuntimeAttributes(job)
    val key = BackendJobDescriptorKey(job, None, 1)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration))
    val testActorRef = TestActorRef[TestablePipelinesApiJobExecutionActor](
      props,
      s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    testActorRef.underlyingActor.monitoringScript shouldBe
      Option(
        PipelinesApiFileInput("monitoring-in",
                              gcsPath("gs://path/to/script"),
                              DefaultPathBuilder.get("monitoring.sh"),
                              workingDisk
        )
      )
  }

  it should "not create a JesFileInput for the monitoring script, when not specified" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow
        .load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
              Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      NoOptions,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val job: CommandCallNode = workflowDescriptor.callable.graph.nodes.collectFirst { case t: CommandCallNode => t }.get
    val key = BackendJobDescriptorKey(job, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(job)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration))
    val testActorRef = TestActorRef[TestablePipelinesApiJobExecutionActor](
      props,
      s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    testActorRef.underlyingActor.monitoringScript shouldBe None
  }

  it should "return JES log paths for non-scattered call" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7c")),
      WdlNamespaceWithWorkflow
        .load(SampleWdl.HelloWorld.asWorkflowSources(""" runtime {docker: "ubuntu:latest"} """).workflowSource.get,
              Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions.fromJsonString(""" {"jes_gcs_root": "gs://path/to/gcs_root"} """).get,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == "hello").get
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration))
    val testActorRef = TestActorRef[TestablePipelinesApiJobExecutionActor](
      props,
      s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    val jesBackend = testActorRef.underlyingActor

    jesBackend.pipelinesApiCallPaths.stdout should be(a[GcsPath])
    jesBackend.pipelinesApiCallPaths.stdout.pathAsString shouldBe
      "gs://path/to/gcs_root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/stdout"
    jesBackend.pipelinesApiCallPaths.stderr should be(a[GcsPath])
    jesBackend.pipelinesApiCallPaths.stderr.pathAsString shouldBe
      "gs://path/to/gcs_root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/stderr"
    jesBackend.pipelinesApiCallPaths.jesLogPath should be(a[GcsPath])
    jesBackend.pipelinesApiCallPaths.jesLogPath.pathAsString shouldBe
      "gs://path/to/gcs_root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello.log"
  }

  it should "return JES log paths for scattered call" in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7d")),
      WdlNamespaceWithWorkflow
        .load(
          new SampleWdl.ScatterWdl().asWorkflowSources(""" runtime {docker: "ubuntu:latest"} """).workflowSource.get,
          Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions.fromJsonString(""" {"jes_gcs_root": "gs://path/to/gcs_root"} """).get,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == "B").get
    val key = BackendJobDescriptorKey(call, Option(2), 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration))
    val testActorRef = TestActorRef[TestablePipelinesApiJobExecutionActor](
      props,
      s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    val jesBackend = testActorRef.underlyingActor

    jesBackend.pipelinesApiCallPaths.stdout should be(a[GcsPath])
    jesBackend.pipelinesApiCallPaths.stdout.pathAsString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/stdout"
    jesBackend.pipelinesApiCallPaths.stderr should be(a[GcsPath])
    jesBackend.pipelinesApiCallPaths.stderr.pathAsString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/stderr"
    jesBackend.pipelinesApiCallPaths.jesLogPath should be(a[GcsPath])
    jesBackend.pipelinesApiCallPaths.jesLogPath.pathAsString shouldBe
      "gs://path/to/gcs_root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/B-2.log"
  }

  it should "return preemptible = true only in the correct cases" in {
    def attempt(max: Int, attempt: Int): PipelinesApiAsyncBackendJobExecutionActor =
      buildPreemptibleTestActorRef(attempt, max).underlyingActor
    def attempt1(max: Int) = attempt(max, 1)
    def attempt2(max: Int) = attempt(max, 2)

    val descriptorWithMax0AndKey1 = attempt1(max = 0)
    descriptorWithMax0AndKey1.preemptible shouldBe false

    val descriptorWithMax1AndKey1 = attempt1(max = 1)
    descriptorWithMax1AndKey1.preemptible shouldBe true

    val descriptorWithMax2AndKey1 = attempt1(max = 2)
    descriptorWithMax2AndKey1.preemptible shouldBe true

    val descriptorWithMax1AndKey2 = attempt2(max = 1)
    descriptorWithMax1AndKey2.preemptible shouldBe false

    val descriptorWithMax2AndKey2 = attempt2(max = 2)
    descriptorWithMax2AndKey2.preemptible shouldBe true
  }

  it should "return the project from the workflow options in the start metadata" in {
    val googleProject = "baa-ram-ewe"
    val jesGcsRoot = "gs://anorexic/duck"
    val workflowId = WorkflowId.randomId()
    val workflowDescriptor = BackendWorkflowDescriptor(
      workflowId,
      WdlNamespaceWithWorkflow
        .load(
          SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
          Seq.empty[Draft2ImportResolver]
        )
        .get
        .workflow
        .toWomWorkflowDefinition(isASubworkflow = false)
        .getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions
        .fromJsonString(
          s"""|{
              |  "google_project": "$googleProject",
              |  "jes_gcs_root": "$jesGcsRoot"
              |}
              |""".stripMargin
        )
        .get,
      Labels.empty,
      HogGroup("foo"),
      List.empty,
      None
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == "goodbye").get
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor =
      BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestablePipelinesApiJobExecutionActor(jobDescriptor, Promise(), papiConfiguration))
    val testActorRef = TestActorRef[TestablePipelinesApiJobExecutionActor](
      props,
      s"TestableJesJobExecutionActor-${jobDescriptor.workflowDescriptor.id}"
    )

    val jesBackend = testActorRef.underlyingActor

    val actual = jesBackend.startMetadataKeyValues.safeMapValues(_.toString)
    actual should be(
      Map(
        "backendLogs:log" -> s"$jesGcsRoot/wf_hello/$workflowId/call-goodbye/goodbye.log",
        "callRoot" -> s"$jesGcsRoot/wf_hello/$workflowId/call-goodbye",
        "jes:endpointUrl" -> "https://lifesciences.googleapis.com/",
        "jes:executionBucket" -> jesGcsRoot,
        "jes:googleProject" -> googleProject,
        "labels:cromwell-workflow-id" -> s"cromwell-$workflowId",
        "labels:wdl-task-name" -> "goodbye",
        "preemptible" -> "false",
        "runtimeAttributes:bootDiskSizeGb" -> "10",
        "runtimeAttributes:continueOnReturnCode" -> "0",
        "runtimeAttributes:cpu" -> "1",
        "runtimeAttributes:disks" -> "local-disk 200 SSD",
        "runtimeAttributes:docker" -> "ubuntu:latest",
        "runtimeAttributes:failOnStderr" -> "false",
        "runtimeAttributes:memory" -> "2 GB",
        "runtimeAttributes:noAddress" -> "false",
        "runtimeAttributes:preemptible" -> "0",
        "runtimeAttributes:zones" -> "us-central1-b,us-central1-a",
        "runtimeAttributes:maxRetries" -> "0",
        "stderr" -> s"$jesGcsRoot/wf_hello/$workflowId/call-goodbye/stderr",
        "stdout" -> s"$jesGcsRoot/wf_hello/$workflowId/call-goodbye/stdout"
      )
    )

  }

  private def makeRuntimeAttributes(job: CommandCallNode) = {
    val evaluatedAttributes =
      RuntimeAttributeDefinition.evaluateRuntimeAttributes(job.callable.runtimeAttributes, null, Map.empty)
    RuntimeAttributeDefinition.addDefaultsToAttributes(runtimeAttributesBuilder.definitions.toSet, NoOptions)(
      evaluatedAttributes.getOrElse(fail("Failed to evaluate runtime attributes"))
    )
  }
}
