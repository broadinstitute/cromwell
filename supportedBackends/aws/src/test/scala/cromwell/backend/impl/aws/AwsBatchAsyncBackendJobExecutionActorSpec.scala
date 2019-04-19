/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

package cromwell.backend.impl.aws

import java.util.UUID

import akka.actor.{ActorRef, Props}
import akka.testkit.{ImplicitSender, TestActorRef, TestDuration}
import common.collections.EnhancedCollections._
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend._
import cromwell.backend.async.{ExecutionHandle, FailedNonRetryableExecutionHandle}
import cromwell.backend.impl.aws.AwsBatchAsyncBackendJobExecutionActor.AwsBatchPendingExecutionHandle
import cromwell.backend.impl.aws.RunStatus.UnsuccessfulRunStatus
import cromwell.backend.impl.aws.io.AwsBatchWorkingDisk
import cromwell.backend.io.JobPathsSpecHelper._
import cromwell.backend.standard.{DefaultStandardAsyncExecutionActorParams, StandardAsyncExecutionActorParams, StandardAsyncJob, StandardExpressionFunctionsParams}
import cromwell.cloudsupport.aws.s3.S3Storage
import cromwell.core.Tags.AwsTest
import cromwell.core._
import cromwell.core.callcaching.NoDocker
import cromwell.core.labels.Labels
import cromwell.core.logging.JobLogger
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import cromwell.filesystems.s3.{S3Path, S3PathBuilder}
import cromwell.services.keyvalue.InMemoryKvServiceActor
import cromwell.services.keyvalue.KeyValueServiceActor.KvResponse
import cromwell.util.JsonFormatting.WomValueJsonFormatter._
import cromwell.util.SampleWdl
import org.scalatest._
import org.slf4j.Logger
import org.specs2.mock.Mockito
import software.amazon.awssdk.auth.credentials.AnonymousCredentialsProvider
import software.amazon.awssdk.regions.Region
import spray.json._
import wdl.draft2.model._
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
import scala.concurrent.{Await, Future, Promise}
import scala.language.postfixOps
import scala.util.Success

class AwsBatchAsyncBackendJobExecutionActorSpec extends TestKitSuite("AwsBatchAsyncBackendJobExecutionActorSpec")
  with FlatSpecLike with Matchers with ImplicitSender with Mockito with BackendSpec with BeforeAndAfter with DefaultJsonProtocol {
  lazy val mockPathBuilder: S3PathBuilder = S3PathBuilder.fromCredentials(
    AnonymousCredentialsProvider.create.resolveCredentials(),
    S3Storage.DefaultConfiguration,
    WorkflowOptions.empty,
    Option(Region.US_EAST_1)
  )

  var kvService: ActorRef = system.actorOf(Props(new InMemoryKvServiceActor))

  import AwsBatchTestConfig._

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
      |    docker: "alpine:latest"
      |    queueArn: "arn:aws:myarn"
      |  }
      |}
      |
      |workflow wf_sup {
      |  call sup
      |}
    """.stripMargin

  val Inputs: Map[FullyQualifiedName, WomValue] = Map("wf_sup.sup.addressee" -> WomString("dog"))

  val NoOptions = WorkflowOptions(JsObject(Map.empty[String, JsValue]))

  lazy val TestableCallContext =
    CallContext(mockPathBuilder.build("s3://root").get, DummyStandardPaths, isDocker = false)

  lazy val TestableStandardExpressionFunctionsParams = new StandardExpressionFunctionsParams {
    override lazy val pathBuilders: List[PathBuilder] = List(mockPathBuilder)
    override lazy val callContext: CallContext = TestableCallContext
    override val ioActorProxy: ActorRef = simpleIoActor
    override val executionContext = system.dispatcher
  }

  lazy val TestableAwsBatchExpressionFunctions: AwsBatchExpressionFunctions = {
    new AwsBatchExpressionFunctions(TestableStandardExpressionFunctionsParams)
  }

  private def buildInitializationData(jobDescriptor: BackendJobDescriptor, configuration: AwsBatchConfiguration) = {
    val workflowPaths = AwsBatchWorkflowPaths(
      jobDescriptor.workflowDescriptor,
      AnonymousCredentialsProvider.create.resolveCredentials(),
      configuration
    )
    val runtimeAttributesBuilder = AwsBatchRuntimeAttributes.runtimeAttributesBuilder(configuration)
    AwsBatchBackendInitializationData(workflowPaths, runtimeAttributesBuilder, configuration, null)
  }

  class TestableAwsBatchJobExecutionActor(params: StandardAsyncExecutionActorParams, functions: AwsBatchExpressionFunctions)
    extends AwsBatchAsyncBackendJobExecutionActor(params) {

    def this(jobDescriptor: BackendJobDescriptor,
             promise: Promise[BackendJobExecutionResponse],
             configuration: AwsBatchConfiguration,
             functions: AwsBatchExpressionFunctions = TestableAwsBatchExpressionFunctions,
             singletonActor: ActorRef = emptyActor,
             ioActor: ActorRef = mockIoActor) = {

      this(
        DefaultStandardAsyncExecutionActorParams(
          jobIdKey = AwsBatchAsyncBackendJobExecutionActor.AwsBatchOperationIdKey,
          serviceRegistryActor = kvService,
          ioActor = ioActor,
          jobDescriptor = jobDescriptor,
          configurationDescriptor = configuration.configurationDescriptor,
          backendInitializationDataOption = Option(buildInitializationData(jobDescriptor, configuration)),
          backendSingletonActorOption = Option(singletonActor),
          completionPromise = promise,
          minimumRuntimeSettings = MinimumRuntimeSettings()
        ),
        functions
      )
    }

    override lazy val jobLogger = new JobLogger(
      "TestLogger",
      workflowIdForLogging,
      rootWorkflowIdForLogging,
      jobTag,
      akkaLogger = Option(log)
    ) {
      override def tag: String = s"$name [UUID(${workflowIdForLogging.shortString})$jobTag]"
      override val slf4jLoggers: Set[Logger] = Set.empty
    }

    override lazy val backendEngineFunctions: AwsBatchExpressionFunctions = functions
  }

  private val configuration = new AwsBatchConfiguration(AwsBatchBackendConfigurationDescriptor)
  private val runtimeAttributesBuilder = AwsBatchRuntimeAttributes.runtimeAttributesBuilder(configuration)
  private val workingDisk = AwsBatchWorkingDisk()

  val DockerAndDiskRuntime: String =
    """
      |runtime {
      |  docker: "alpine:latest"
      |  disks: "local-disk"
      |  queueArn: "arn:aws:myarn"
      |}
    """.stripMargin

  private def buildJobDescriptor(): BackendJobDescriptor = {
    val attempt = 1
    val wdlNamespace = WdlNamespaceWithWorkflow.load(YoSup, Seq.empty[Draft2ImportResolver]).get
    val womDefinition = wdlNamespace.workflow.toWomWorkflowDefinition(isASubworkflow = false).getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))

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
          Labels.empty
        )

        val job = workflowDescriptor.callable.taskCallNodes.head
        val key = BackendJobDescriptorKey(job, None, attempt)
        val runtimeAttributes = makeRuntimeAttributes(job)
        val prefetchedKvEntries = Map[String, KvResponse]()
        BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(Inputs), NoDocker, None, prefetchedKvEntries)
      case Left(badtimes) => fail(badtimes.toList.mkString(", "))
    }
  }

  // private def executionActor(jobDescriptor: BackendJobDescriptor,
  //                            configurationDescriptor: BackendConfigurationDescriptor,
  //                            promise: Promise[BackendJobExecutionResponse],
  //                            singletonActor: ActorRef
  //                            ): ActorRef = {
  //
  //   val job = StandardAsyncJob(UUID.randomUUID().toString)
  //   val run = Run(job)
  //   val handle = new AwsBatchPendingExecutionHandle(jobDescriptor, run.job, Option(run), None)
  //
  //   class ExecuteOrRecoverActor extends TestableAwsBatchJobExecutionActor(jobDescriptor, promise, configuration, singletonActor = singletonActor) {
  //     override def executeOrRecover(mode: ExecutionMode)(implicit ec: ExecutionContext): Future[ExecutionHandle] = {
  //       Future.successful(handle)
  //     }
  //   }
  //
  //   system.actorOf(Props(new ExecuteOrRecoverActor), "ExecuteOrRecoverActor-" + UUID.randomUUID)
  // }

  // private def runAndFail(previousPreemptions: Int, previousUnexpectedRetries: Int, errorCode: Status, innerErrorMessage: String): BackendJobExecutionResponse = {
  //
  //   // val runStatus = UnsuccessfulRunStatus("test", "failed", errorCode, Option(innerErrorMessage), Seq.empty)
  //   val statusPoller = TestProbe()
  //
  //   val promise = Promise[BackendJobExecutionResponse]()
  //   val jobDescriptor =  buildJobDescriptor()
  //
  //   // TODO: Use this to check the new KV entries are there!
  //   //val kvProbe = TestProbe()
  //
  //   val backend = executionActor(jobDescriptor, AwsBatchBackendConfigurationDescriptor, promise, statusPoller.ref)
  //   backend ! Execute
  //   // statusPoller.expectMsgPF(max = Timeout, hint = "awaiting status poll") {
  //   //   case _: DoPoll => backend ! runStatus
  //   // }
  //
  //   Await.result(promise.future, Timeout)
  // }

  def buildTestActorRef(attempt: Int): TestActorRef[TestableAwsBatchJobExecutionActor] = {
    // For this test we say that all previous attempts were preempted:
    val jobDescriptor = buildJobDescriptor()
    val props = Props(new TestableAwsBatchJobExecutionActor(jobDescriptor, Promise(),
      configuration,
      TestableAwsBatchExpressionFunctions,
      emptyActor,
      failIoActor))
    TestActorRef(props, s"TestableAwsBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")
  }

  behavior of "AwsBatchAsyncBackendJobExecutionActor"

  val timeout = 25 seconds

//   { // Set of "handle call failures appropriately with respect to preemption and failure" tests
//     val expectations = Table(
//       ("previous_preemptions", "previous_unexpectedRetries", "errorCode", "message", "shouldRetry"),
//       // No preemptible attempts allowed, but standard failures should be retried.
//       (0, 0, Status.ABORTED, "13: retryable error", true), // This is the new "unexpected failure" mode, which is now retried
//       (0, 1, Status.ABORTED, "13: retryable error", true),
//       (0, 2, Status.ABORTED, "13: retryable error", false), // The third unexpected failure is a real failure.
//       (0, 0, Status.ABORTED, "14: usually means preempted...?", false), // Usually means "preempted', but this wasn't a preemptible VM, so this should just be a failure.
//       (0, 0, Status.ABORTED, "15: other error", false),
//       (0, 0, Status.OUT_OF_RANGE, "13: unexpected error", false),
//       (0, 0, Status.OUT_OF_RANGE, "14: test error msg", false),
//       // These commented out tests should be uncommented if/when we stop mapping 13 to 14 in preemption mode
//       // 1 preemptible attempt allowed, but not all failures represent preemptions.
// //      (0, 0, Status.ABORTED, "13: retryable error", true),
// //      (0, 1, Status.ABORTED, "13: retryable error", true),
// //      (0, 2, Status.ABORTED, "13: retryable error", false),
//       // The following 13 based test should be removed if/when we stop mapping 13 to 14 in preemption mode
//       (0, 0, Status.ABORTED, "14: preempted", true),
//       (0, 0, Status.UNKNOWN, "Instance failed to start due to preemption.", true),
//       (0, 0, Status.OUT_OF_RANGE, "13: retryable error", false),
//       (0, 0, Status.OUT_OF_RANGE, "14: preempted", false),
//       (0, 0, Status.OUT_OF_RANGE, "Instance failed to start due to preemption.", false),
//       // 1 preemptible attempt allowed, but since we're now on the second preemption attempt only 13s should be retryable.
//       (1, 0, Status.ABORTED, "13: retryable error", true),
//       (1, 1, Status.ABORTED, "13: retryable error", true),
//       (1, 2, Status.ABORTED, "13: retryable error", false),
//       (1, 0, Status.ABORTED, "14: preempted", false),
//       (1, 0, Status.UNKNOWN, "Instance failed to start due to preemption.", false),
//       (1, 0, Status.ABORTED, "15: other error", false),
//       (1, 0, Status.OUT_OF_RANGE, "13: retryable error", false),
//       (1, 0, Status.OUT_OF_RANGE, "14: preempted", false),
//       (1, 0, Status.OUT_OF_RANGE, "Instance failed to start due to preemption.", false)
//     )
//
//     expectations foreach { case (previousPreemptions, previousUnexpectedRetries, errorCode, innerErrorMessage, shouldRetry) =>
//       val descriptor = s"previousPreemptions=$previousPreemptions, previousUnexpectedRetries=$previousUnexpectedRetries errorCode=$errorCode, innerErrorMessage=$innerErrorMessage"
//       it should s"handle call failures appropriately with respect to preemption and failure ($descriptor)" in {
//         runAndFail(previousPreemptions, previousUnexpectedRetries, errorCode, innerErrorMessage) match {
//           case response: JobFailedNonRetryableResponse =>
//             if(shouldRetry) fail(s"A should-be-retried job ($descriptor) was sent back to the engine with: $response")
//           case response: JobFailedRetryableResponse =>
//             if(!shouldRetry) fail(s"A shouldn't-be-retried job ($descriptor) was sent back to the engine with $response")
//           case huh => fail(s"Unexpected response: $huh")
//         }
//       }
//     }
//   }

  // TODO: Rework all this stuff
  // it should "not restart 2 of 1 unexpected shutdowns without another preemptible VM" in {
  //   val actorRef = buildTestActorRef(1)
  //   val backend = actorRef.underlyingActor
  //   val runId = StandardAsyncJob(UUID.randomUUID().toString)
  //   val handle = new AwsBatchPendingExecutionHandle(null, runId, None, None)
  //
  //   // TODO: Get real error messages
  //   val failedStatus = UnsuccessfulRunStatus("test", "200", Status.ABORTED, Option("14: VM XXX shut down unexpectedly."), Seq.empty)
  //   val executionResult = backend.handleExecutionResult(failedStatus, handle)
  //   val result = Await.result(executionResult, timeout)
  //   result.isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
  //   val failedHandle = result.asInstanceOf[FailedNonRetryableExecutionHandle]
  //   failedHandle.returnCode shouldBe None
  // }
  //
  // it should "restart 1 of 1 unexpected shutdowns without another preemptible VM" in {
  //   val actorRef = buildTestActorRef(1)
  //   val backend = actorRef.underlyingActor
  //   val runId = StandardAsyncJob(UUID.randomUUID().toString)
  //   val handle = new AwsBatchPendingExecutionHandle(null, runId, None, None)
  //
  //   val failedStatus = UnsuccessfulRunStatus("test", "200", Status.ABORTED, Option("14: VM XXX shut down unexpectedly."), Seq.empty)
  //   val executionResult = backend.handleExecutionResult(failedStatus, handle)
  //   val result = Await.result(executionResult, timeout)
  //   result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
  //   val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
  //   retryableHandle.returnCode shouldBe None
  //   retryableHandle.throwable.getMessage should include("will be restarted with a non-preemptible VM")
  // }
  //
  // it should "restart 1 of 2 unexpected shutdowns with another preemptible VM" in {
  //   val actorRef = buildTestActorRef(2)
  //   val backend = actorRef.underlyingActor
  //   val runId = StandardAsyncJob(UUID.randomUUID().toString)
  //   val handle = new AwsBatchPendingExecutionHandle(null, runId, None, None)
  //
  //   val failedStatus = UnsuccessfulRunStatus("test", "200", Status.ABORTED, Option("14: VM XXX shut down unexpectedly."), Seq.empty)
  //   val executionResult = backend.handleExecutionResult(failedStatus, handle)
  //   val result = Await.result(executionResult, timeout)
  //   result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
  //   val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
  //   retryableHandle.returnCode shouldBe None
  //   retryableHandle.throwable.getMessage should include("will be restarted with another preemptible VM")
  // }
  //
  // it should "handle message 13 (TODO: fix this)" in {
  //   val actorRef = buildTestActorRef(2)
  //   val backend = actorRef.underlyingActor
  //   val runId = StandardAsyncJob(UUID.randomUUID().toString)
  //   val handle = new AwsBatchPendingExecutionHandle(null, runId, None, None)
  //
  //   val failedStatus = UnsuccessfulRunStatus("test", "200", Status.ABORTED, Option("13: Retryable Error."), Seq.empty)
  //   val executionResult = backend.handleExecutionResult(failedStatus, handle)
  //   val result = Await.result(executionResult, timeout)
  //   result.isInstanceOf[FailedRetryableExecutionHandle] shouldBe true
  //   val retryableHandle = result.asInstanceOf[FailedRetryableExecutionHandle]
  //   retryableHandle.returnCode shouldBe None
  //   retryableHandle.throwable.getMessage should include("will be restarted with another preemptible VM")
  // }

  it should "handle Failure Status for various errors" taggedAs AwsTest in {
    val actorRef = buildTestActorRef(1)
    val backend = actorRef.underlyingActor
    val runId = StandardAsyncJob(UUID.randomUUID().toString)
    val handle = new AwsBatchPendingExecutionHandle(null, runId, None, None)

    def checkFailedResult(errorMessage: Option[String]): ExecutionHandle = {
      val failed = UnsuccessfulRunStatus("test", "failed", errorMessage, Seq.empty)
      Await.result(backend.handleExecutionResult(failed, handle), timeout)
    }

    checkFailedResult(Option("15: Other type of error."))
      .isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Option("14: Wrong errorCode.")).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Option("Weird error message.")).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(Option("UnparsableInt: Even weirder error message."))
      .isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    checkFailedResult(None).isInstanceOf[FailedNonRetryableExecutionHandle] shouldBe true
    // TODO: Determine the actual error message that comes back and special-case it in handlExecutionFailure
    //       in AwsBatchBackendJobExecutionActor
    //checkFailedResult(Option("Operation canceled at")) shouldBe AbortedExecutionHandle

    actorRef.stop()
  }

  it should "map paths and *only* foreign paths to local" taggedAs AwsTest ignore {
    val stringKey = "abc"
    val stringVal = WomString("abc")
    val localFileKey = "lf"
    val localFileVal = WomSingleFile("/blah/abc")
    val fileKey = "batchf"
    val fileVal = WomSingleFile("s3://blah/abc")

    val inputs: Map[String, WomValue] = Map(
      stringKey -> stringVal,
      localFileKey -> localFileVal,
      fileKey -> fileVal
    )

    val wdlNamespace = WdlNamespaceWithWorkflow.load(YoSup,
      Seq.empty[Draft2ImportResolver]).get
    val womWorkflow = wdlNamespace.workflow.toWomWorkflowDefinition(isASubworkflow = false).getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    wdlNamespace.toWomExecutable(Option(inputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }})

        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty
        )

        val call: CommandCallNode = workflowDescriptor.callable.graph.nodes.collectFirst({ case t: CommandCallNode => t }).get
        val key = BackendJobDescriptorKey(call, None, 1)
        val runtimeAttributes = makeRuntimeAttributes(call)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, None, Map.empty)

        val props = Props(new TestableAwsBatchJobExecutionActor(jobDescriptor, Promise(), configuration))
        val testActorRef = TestActorRef[TestableAwsBatchJobExecutionActor](
          props, s"TestableAwsBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")


        def pathToLocal(womValue: WomValue): WomValue = {
          WomFileMapper.mapWomFiles(testActorRef.underlyingActor.mapCommandLineWomFile, exceptions = Set.empty)(womValue).get
        }

        val mappedInputs = jobDescriptor.localInputs safeMapValues pathToLocal

        mappedInputs(stringKey) match {
          case WomString(v) => assert(v.equalsIgnoreCase(stringVal.value))
          case _ => fail("test setup error")
        }

        mappedInputs(localFileKey) match {
          case wdlFile: WomSingleFile => assert(wdlFile.value.equalsIgnoreCase(localFileVal.value))
          case _ => fail("test setup error")
        }

        mappedInputs(fileKey) match {
          case wdlFile: WomSingleFile => assert(wdlFile.value.equalsIgnoreCase("/cromwell_root/blah/abc"))
          case _ => fail("test setup error")
        }
      case Left(badtimes) => fail(badtimes.toList.mkString(", "))
    }
  }

  private val dockerAndDiskWdlNamespace = WdlNamespaceWithWorkflow.load(SampleWdl.CurrentDirectory.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
    Seq.empty[Draft2ImportResolver]).get

  it should "generate correct AwsBatchFileInputs from a WdlMap" taggedAs AwsTest ignore {
    val inputs: Map[String, WomValue] = Map(
      "stringToFileMap" -> WomMap(WomMapType(WomStringType, WomSingleFileType), Map(
        WomString("stringTofile1") -> WomSingleFile("s3://path/to/stringTofile1"),
        WomString("stringTofile2") -> WomSingleFile("s3://path/to/stringTofile2")
      )),
      "fileToStringMap" -> WomMap(WomMapType(WomSingleFileType, WomStringType), Map(
        WomSingleFile("s3://path/to/fileToString1") -> WomString("fileToString1"),
        WomSingleFile("s3://path/to/fileToString2") -> WomString("fileToString2")
      )),
      "fileToFileMap" -> WomMap(WomMapType(WomSingleFileType, WomSingleFileType), Map(
        WomSingleFile("s3://path/to/fileToFile1Key") -> WomSingleFile("s3://path/to/fileToFile1Value"),
        WomSingleFile("s3://path/to/fileToFile2Key") -> WomSingleFile("s3://path/to/fileToFile2Value")
      )),
      "stringToString" -> WomMap(WomMapType(WomStringType, WomStringType), Map(
        WomString("stringToString1") -> WomString("path/to/stringToString1"),
        WomString("stringToString2") -> WomString("path/to/stringToString2")
      ))
    )

    val workflowInputs = inputs map {
      case (k, v) => s"wf_whereami.whereami$k" -> v
    }

    val womWorkflow = dockerAndDiskWdlNamespace.workflow.toWomWorkflowDefinition(isASubworkflow = false).getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    dockerAndDiskWdlNamespace.toWomExecutable(Option(workflowInputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }})
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, None, Map.empty)

        val props = Props(new TestableAwsBatchJobExecutionActor(jobDescriptor, Promise(), configuration))
        val testActorRef = TestActorRef[TestableAwsBatchJobExecutionActor](
          props, s"TestableAwsBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

        val batchInputs = testActorRef.underlyingActor.generateAwsBatchInputs(jobDescriptor)
        batchInputs should have size 8
        batchInputs should contain(AwsBatchFileInput(
          "stringToFileMap-0", "s3://path/to/stringTofile1", DefaultPathBuilder.get("path/to/stringTofile1"), workingDisk))
        batchInputs should contain(AwsBatchFileInput(
          "stringToFileMap-1", "s3://path/to/stringTofile2", DefaultPathBuilder.get("path/to/stringTofile2"), workingDisk))
        batchInputs should contain(AwsBatchFileInput(
          "fileToStringMap-0", "s3://path/to/fileToString1", DefaultPathBuilder.get("path/to/fileToString1"), workingDisk))
        batchInputs should contain(AwsBatchFileInput(
          "fileToStringMap-1", "s3://path/to/fileToString2", DefaultPathBuilder.get("path/to/fileToString2"), workingDisk))
        batchInputs should contain(AwsBatchFileInput(
          "fileToFileMap-0", "s3://path/to/fileToFile1Key", DefaultPathBuilder.get("path/to/fileToFile1Key"), workingDisk))
        batchInputs should contain(AwsBatchFileInput(
          "fileToFileMap-1", "s3://path/to/fileToFile1Value", DefaultPathBuilder.get("path/to/fileToFile1Value"), workingDisk))
        batchInputs should contain(AwsBatchFileInput(
          "fileToFileMap-2", "s3://path/to/fileToFile2Key", DefaultPathBuilder.get("path/to/fileToFile2Key"), workingDisk))
        batchInputs should contain(AwsBatchFileInput(
          "fileToFileMap-3", "s3://path/to/fileToFile2Value", DefaultPathBuilder.get("path/to/fileToFile2Value"), workingDisk))

      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  def makeAwsBatchActorRef(sampleWdl: SampleWdl, callName: LocallyQualifiedName, inputs: Map[FullyQualifiedName, WomValue],
                      functions: AwsBatchExpressionFunctions = TestableAwsBatchExpressionFunctions):
  TestActorRef[TestableAwsBatchJobExecutionActor] = {
    val womWorkflow = WdlNamespaceWithWorkflow.load(sampleWdl.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
      Seq.empty[Draft2ImportResolver]).get.workflow.toWomWorkflowDefinition(isASubworkflow = false).getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    dockerAndDiskWdlNamespace.toWomExecutable(Option(inputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }})
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty
        )

        val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == callName).get
        val key = BackendJobDescriptorKey(call, None, 1)
        val runtimeAttributes = makeRuntimeAttributes(call)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, None, Map.empty)

        val props = Props(new TestableAwsBatchJobExecutionActor(jobDescriptor, Promise(), configuration, functions))
        TestActorRef[TestableAwsBatchJobExecutionActor](props, s"TestableAwsBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")
      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "generate correct AwsBatchOutputs" taggedAs AwsTest ignore {
    val inputs = Map(
      "in" -> WomSingleFile("s3://blah/b/c.txt")
    )
    val backend = makeAwsBatchActorRef(SampleWdl.FilePassingWorkflow, "a", inputs).underlyingActor
    val jobDescriptor = backend.jobDescriptor
    val workflowId = backend.workflowDescriptor.id
    val batchInputs = backend.generateAwsBatchInputs(jobDescriptor)
    batchInputs should have size 1
    batchInputs should contain(AwsBatchFileInput("in-0", "s3://blah/b/c.txt", DefaultPathBuilder.get("blah/b/c.txt"), workingDisk))
    val outputs = backend.generateAwsBatchOutputs(jobDescriptor)
    outputs should have size 1
    outputs should contain(AwsBatchFileOutput("out",
      s"s3://my-cromwell-workflows-bucket/file_passing/$workflowId/call-a/out", DefaultPathBuilder.get("out"), workingDisk))
  }

  it should "generate correct AwsBatchInputs when a command line contains a write_lines call in it" taggedAs AwsTest ignore {
    val inputs = Map(
      "strs" -> WomArray(WomArrayType(WomStringType), Seq("A", "B", "C").map(WomString))
    )

    class TestAwsBatchExpressionFunctions extends AwsBatchExpressionFunctions(TestableStandardExpressionFunctionsParams) {
      override def writeFile(path: String, content: String): Future[WomSingleFile] = {
        Future.fromTry(Success(WomSingleFile(s"s3://some/path/file.txt")))
      }
    }

    val functions = new TestAwsBatchExpressionFunctions
    val backend = makeAwsBatchActorRef(SampleWdl.ArrayIO, "serialize", inputs, functions).underlyingActor
    val jobDescriptor = backend.jobDescriptor
    val batchInputs = backend.generateAwsBatchInputs(jobDescriptor)
    batchInputs should have size 1
    batchInputs should contain(AwsBatchFileInput(
      "c6fd5c91-0", "s3://some/path/file.txt", DefaultPathBuilder.get("some/path/file.txt"), workingDisk))
    val outputs = backend.generateAwsBatchOutputs(jobDescriptor)
    outputs should have size 0
  }

  it should "generate correct AwsBatchFileInputs from a WdlArray" taggedAs AwsTest ignore {
    val inputs: Map[String, WomValue] = Map(
      "fileArray" ->
        WomArray(WomArrayType(WomSingleFileType), Seq(WomSingleFile("s3://path/to/file1"), WomSingleFile("s3://path/to/file2")))
    )

    val womWorkflow = dockerAndDiskWdlNamespace.workflow.toWomWorkflowDefinition(isASubworkflow = false).getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    dockerAndDiskWdlNamespace.toWomExecutable(Option(inputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }})
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, None, Map.empty)

        val props = Props(new TestableAwsBatchJobExecutionActor(jobDescriptor, Promise(), configuration))
        val testActorRef = TestActorRef[TestableAwsBatchJobExecutionActor](
          props, s"TestableAwsBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

        val batchInputs = testActorRef.underlyingActor.generateAwsBatchInputs(jobDescriptor)
        batchInputs should have size 2
        batchInputs should contain(AwsBatchFileInput("fileArray-0", "s3://path/to/file1", DefaultPathBuilder.get("path/to/file1"), workingDisk))
        batchInputs should contain(AwsBatchFileInput("fileArray-1", "s3://path/to/file2", DefaultPathBuilder.get("path/to/file2"), workingDisk))
      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "generate correct AwsBatchFileInputs from a WdlFile" taggedAs AwsTest ignore {
    val inputs: Map[String, WomValue] = Map(
      "file1" -> WomSingleFile("s3://path/to/file1"),
      "file2" -> WomSingleFile("s3://path/to/file2")
    )

    val womWorkflow = dockerAndDiskWdlNamespace.workflow.toWomWorkflowDefinition(isASubworkflow = false).getOrElse(fail("failed to get WomDefinition from WdlWorkflow"))
    dockerAndDiskWdlNamespace.toWomExecutable(Option(inputs.toJson.compactPrint), NoIoFunctionSet, strictValidation = true) match {
      case Right(womExecutable) =>
        val wdlInputs = womExecutable.resolvedExecutableInputs.flatMap({case (port, v) => v.select[WomValue] map { port -> _ }})
        val workflowDescriptor = BackendWorkflowDescriptor(
          WorkflowId.randomId(),
          womWorkflow,
          wdlInputs,
          NoOptions,
          Labels.empty
        )

        val job: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
        val runtimeAttributes = makeRuntimeAttributes(job)
        val key = BackendJobDescriptorKey(job, None, 1)
        val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, fqnWdlMapToDeclarationMap(inputs), NoDocker, None, Map.empty)

        val props = Props(new TestableAwsBatchJobExecutionActor(jobDescriptor, Promise(), configuration))
        val testActorRef = TestActorRef[TestableAwsBatchJobExecutionActor](
          props, s"TestableAwsBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

        val batchInputs = testActorRef.underlyingActor.generateAwsBatchInputs(jobDescriptor)
        batchInputs should have size 2
        batchInputs should contain(AwsBatchFileInput("file1-0", "s3://path/to/file1", DefaultPathBuilder.get("path/to/file1"), workingDisk))
        batchInputs should contain(AwsBatchFileInput("file2-0", "s3://path/to/file2", DefaultPathBuilder.get("path/to/file2"), workingDisk))

      case Left(badness) => fail(badness.toList.mkString(", "))
    }
  }

  it should "convert local Paths back to corresponding foreign paths in AwsBatchOutputs" taggedAs AwsTest in {
    val outputs = Set(
      AwsBatchFileOutput("/cromwell_root/path/to/file1", "s3://path/to/file1",
        DefaultPathBuilder.get("/cromwell_root/path/to/file1"), workingDisk),
      AwsBatchFileOutput("/cromwell_root/path/to/file2", "s3://path/to/file2",
        DefaultPathBuilder.get("/cromwell_root/path/to/file2"), workingDisk),
      AwsBatchFileOutput("/cromwell_root/path/to/file3", "s3://path/to/file3",
        DefaultPathBuilder.get("/cromwell_root/path/to/file3"), workingDisk),
      AwsBatchFileOutput("/cromwell_root/path/to/file4", "s3://path/to/file4",
        DefaultPathBuilder.get("/cromwell_root/path/to/file4"), workingDisk),
      AwsBatchFileOutput("/cromwell_root/path/to/file5", "s3://path/to/file5",
        DefaultPathBuilder.get("/cromwell_root/path/to/file5"), workingDisk)
    )
    val outputValues = Seq(
      WomSingleFile("/cromwell_root/path/to/file1"),
      WomArray(WomArrayType(WomSingleFileType), Seq(
        WomSingleFile("/cromwell_root/path/to/file2"), WomSingleFile("/cromwell_root/path/to/file3"))),
      WomMap(WomMapType(WomSingleFileType, WomSingleFileType), Map(
        WomSingleFile("/cromwell_root/path/to/file4") -> WomSingleFile("/cromwell_root/path/to/file5")
      ))
    )

    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId.randomId(),
      WdlNamespaceWithWorkflow.load(SampleWdl.EmptyString.asWorkflowSources(DockerAndDiskRuntime).workflowSource.get,
        Seq.empty[Draft2ImportResolver]).get.workflow.toWomWorkflowDefinition(isASubworkflow = false).getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      NoOptions,
      Labels.empty
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.head
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestableAwsBatchJobExecutionActor(jobDescriptor, Promise(), configuration))
    val testActorRef = TestActorRef[TestableAwsBatchJobExecutionActor](
      props, s"TestableAwsBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    def wdlValueToS3Path(outputs: Set[AwsBatchFileOutput])(womValue: WomValue): WomValue = {
      WomFileMapper.mapWomFiles(testActorRef.underlyingActor.womFileToPath(outputs), exceptions = Set.empty)(womValue).get
    }

    val result = outputValues map wdlValueToS3Path(outputs)
    result should have size 3
    result should contain(WomSingleFile("s3://path/to/file1"))
    result should contain(WomArray(WomArrayType(WomSingleFileType),
      Seq(WomSingleFile("s3://path/to/file2"), WomSingleFile("s3://path/to/file3")))
    )
    result should contain(WomMap(WomMapType(WomSingleFileType, WomSingleFileType),
      Map(WomSingleFile("s3://path/to/file4") -> WomSingleFile("s3://path/to/file5")))
    )
  }

  it should "return log paths for non-scattered call" taggedAs AwsTest in {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7c")),
      WdlNamespaceWithWorkflow.load(
        SampleWdl.HelloWorld.asWorkflowSources(""" runtime {docker: "ubuntu:latest"} """).workflowSource.get,
        Seq.empty[Draft2ImportResolver]).get.workflow.toWomWorkflowDefinition(isASubworkflow = false).getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions.fromJsonString(""" {"aws_s3_root": "s3://path/to/root"} """).get,
      Labels.empty
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == "hello").get
    val key = BackendJobDescriptorKey(call, None, 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestableAwsBatchJobExecutionActor(jobDescriptor, Promise(), configuration))
    val testActorRef = TestActorRef[TestableAwsBatchJobExecutionActor](
      props, s"TestableAwsBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val backend = testActorRef.underlyingActor

    backend.callPaths.stdout should be(a[S3Path])
    backend.callPaths.stdout.pathAsString shouldBe
      "s3://path/to/root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello-stdout.log"
    backend.callPaths.stderr should be(a[S3Path])
    backend.callPaths.stderr.pathAsString shouldBe
      "s3://path/to/root/wf_hello/e6236763-c518-41d0-9688-432549a8bf7c/call-hello/hello-stderr.log"
  }

  it should "return log paths for scattered call" taggedAs AwsTest ignore {
    val workflowDescriptor = BackendWorkflowDescriptor(
      WorkflowId(UUID.fromString("e6236763-c518-41d0-9688-432549a8bf7d")),
      WdlNamespaceWithWorkflow.load(
        new SampleWdl.ScatterWdl().asWorkflowSources(""" runtime {docker: "ubuntu:latest"} """).workflowSource.get,
        Seq.empty[Draft2ImportResolver]).get.workflow.toWomWorkflowDefinition(isASubworkflow = false).getOrElse(fail("failed to get WomDefinition from WdlWorkflow")),
      Map.empty,
      WorkflowOptions.fromJsonString(""" {"root": "s3://path/to/root"} """).get,
      Labels.empty
    )

    val call: CommandCallNode = workflowDescriptor.callable.taskCallNodes.find(_.localName == "B").get
    val key = BackendJobDescriptorKey(call, Option(2), 1)
    val runtimeAttributes = makeRuntimeAttributes(call)
    val jobDescriptor = BackendJobDescriptor(workflowDescriptor, key, runtimeAttributes, Map.empty, NoDocker, None, Map.empty)

    val props = Props(new TestableAwsBatchJobExecutionActor(jobDescriptor, Promise(), configuration))
    val testActorRef = TestActorRef[TestableAwsBatchJobExecutionActor](
      props, s"TestableAwsBatchJobExecutionActor-${jobDescriptor.workflowDescriptor.id}")

    val backend = testActorRef.underlyingActor

    backend.callPaths.stdout should be(a[S3Path])
    backend.callPaths.stdout.pathAsString shouldBe
      "s3://path/to/root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/B-2-stdout.log"
    backend.callPaths.stderr should be(a[S3Path])
    backend.callPaths.stderr.pathAsString shouldBe
      "s3://path/to/root/w/e6236763-c518-41d0-9688-432549a8bf7d/call-B/shard-2/B-2-stderr.log"
  }

  private def makeRuntimeAttributes(job: CommandCallNode) = {
    val evaluatedAttributes = RuntimeAttributeDefinition.evaluateRuntimeAttributes(job.callable.runtimeAttributes, null, Map.empty)
    RuntimeAttributeDefinition.addDefaultsToAttributes(
      runtimeAttributesBuilder.definitions.toSet, NoOptions)(evaluatedAttributes.getOrElse(fail("Failed to evaluate runtime attributes")))
  }
}
