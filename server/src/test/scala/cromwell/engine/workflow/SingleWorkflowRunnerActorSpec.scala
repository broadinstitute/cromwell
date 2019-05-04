package cromwell.engine.workflow

import java.time.OffsetDateTime

import akka.actor._
import akka.pattern.ask
import akka.stream.ActorMaterializer
import akka.testkit._
import akka.util.Timeout
import cromwell.CromwellTestKitSpec._
import cromwell._
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.core.{SimpleIoActor, WorkflowSourceFilesCollection}
import cromwell.engine.MockCromwellTerminator
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.SingleWorkflowRunnerActor.RunWorkflow
import cromwell.engine.workflow.SingleWorkflowRunnerActorSpec._
import cromwell.engine.workflow.tokens.DynamicRateLimiter.Rate
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor
import cromwell.engine.workflow.workflowstore._
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.{ExpressionsInInputs, GoodbyeWorld, ThreeStep}
import mouse.all._
import org.scalatest.prop.{TableDrivenPropertyChecks, TableFor3}
import org.specs2.mock.Mockito
import spray.json._

import scala.concurrent.Await
import scala.concurrent.duration._
import scala.util._
import scala.util.control.NoStackTrace

/**
 * A series of tests of the SingleWorkflowRunnerActor. Currently uses live versions of the SingleWorkflowRunnerActor and
 * the WorkflowManagerActor communicating with each other, instead of TestActor/TestProbe.
 *
 * Currently, as instance of the actor system are created via an instance of CromwellTestkitSpec, and the
 * SingleWorkflowRunnerActor also tests halting its actor system, each spec is currently in a separate instance of the
 * CromwellTestKitSpec.
 */
object SingleWorkflowRunnerActorSpec {

  def tempFile() = DefaultPathBuilder.createTempFile("metadata.", ".json")

  def tempDir() = DefaultPathBuilder.createTempDirectory("metadata.dir.")

  implicit class OptionJsValueEnhancer(val jsValue: Option[JsValue]) extends AnyVal {
    def toOffsetDateTime = OffsetDateTime.parse(jsValue.toStringValue)
    def toStringValue = jsValue.getOrElse(JsString("{}")).asInstanceOf[JsString].value
    def toFields = jsValue.get.asJsObject.fields
  }

  class TestSingleWorkflowRunnerActor(source: WorkflowSourceFilesCollection,
                                      metadataOutputPath: Option[Path])(implicit materializer: ActorMaterializer)
    extends SingleWorkflowRunnerActor(
      source = source,
      metadataOutputPath = metadataOutputPath,
      terminator = MockCromwellTerminator,
      gracefulShutdown = false,
      abortJobsOnTerminate = false,
      config = CromwellTestKitSpec.DefaultConfig
    ) {
    override lazy val serviceRegistryActor = CromwellTestKitSpec.ServiceRegistryActorInstance
    override private [workflow] def done() = context.stop(self)
  }
}

abstract class SingleWorkflowRunnerActorSpec extends CromwellTestKitWordSpec with CoordinatedWorkflowStoreBuilder with Mockito {
  private val workflowHeartbeatConfig = WorkflowHeartbeatConfig(CromwellTestKitSpec.DefaultConfig)
  val store = new InMemoryWorkflowStore
  private val workflowStore =
    system.actorOf(
      WorkflowStoreActor.props(
        store,
        store |> access,
        dummyServiceRegistryActor,
        MockCromwellTerminator,
        abortAllJobsOnTerminate = false,
        workflowHeartbeatConfig
      ),
      "WorkflowStoreActor"
    )
  private val serviceRegistry = TestProbe("ServiceRegistryProbe").ref
  private val jobStore = system.actorOf(AlwaysHappyJobStoreActor.props, "AlwaysHappyJobStoreActor")
  private val ioActor = system.actorOf(SimpleIoActor.props, "SimpleIoActor")
  private val subWorkflowStore = system.actorOf(
    AlwaysHappySubWorkflowStoreActor.props,
    "AlwaysHappySubWorkflowStoreActor"
  )
  private val callCacheReadActor = system.actorOf(EmptyCallCacheReadActor.props, "EmptyCallCacheReadActor")
  private val callCacheWriteActor = system.actorOf(EmptyCallCacheWriteActor.props, "EmptyCallCacheWriteActor")
  private val dockerHashActor = system.actorOf(EmptyDockerHashActor.props, "EmptyDockerHashActor")
  private val jobTokenDispenserActor = system.actorOf(
    JobExecutionTokenDispenserActor.props(serviceRegistry, Rate(100, 1.second), None),
    "JobExecutionTokenDispenserActor"
  )


  def workflowManagerActor(): ActorRef = {
    val params = WorkflowManagerActorParams(
      CromwellTestKitSpec.DefaultConfig,
      workflowStore = workflowStore,
      ioActor = ioActor,
      serviceRegistryActor = dummyServiceRegistryActor,
      workflowLogCopyRouter = dummyLogCopyRouter,
      jobStoreActor = jobStore,
      subWorkflowStoreActor = subWorkflowStore,
      callCacheReadActor = callCacheReadActor,
      callCacheWriteActor = callCacheWriteActor,
      dockerHashActor = dockerHashActor,
      jobTokenDispenserActor = jobTokenDispenserActor,
      backendSingletonCollection = BackendSingletonCollection(Map.empty),
      serverMode = false,
      workflowHeartbeatConfig)
    system.actorOf(Props(new WorkflowManagerActor(params)), "WorkflowManagerActor")
  }
  
  def createRunnerActor(sampleWdl: SampleWdl = ThreeStep, managerActor: => ActorRef = workflowManagerActor(),
                          outputFile: => Option[Path] = None): ActorRef = {
    system.actorOf(
      Props(new TestSingleWorkflowRunnerActor(sampleWdl.asWorkflowSources(), outputFile)),
      "TestSingleWorkflowRunnerActor"
    )
  }

  def singleWorkflowActor(sampleWdl: SampleWdl = ThreeStep, managerActor: => ActorRef = workflowManagerActor(),
                          outputFile: => Option[Path] = None): Unit = {
    val actorRef = createRunnerActor(sampleWdl, managerActor, outputFile)
    val futureResult = actorRef.ask(RunWorkflow)(timeout = new Timeout(TimeoutDuration))
    Await.ready(futureResult, Duration.Inf)
    ()
  }
}

class SingleWorkflowRunnerActorNormalSpec extends SingleWorkflowRunnerActorSpec {
  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow" in {
      within(TimeoutDuration) {
        waitForInfo("workflow finished with status 'Succeeded'.") {
          singleWorkflowActor()
        }
      }
    }
  }
}

class SingleWorkflowRunnerActorWithMetadataSpec extends SingleWorkflowRunnerActorSpec with TableDrivenPropertyChecks {
  val metadataFile = tempFile()

  override protected def afterAll() = {
    metadataFile.delete(swallowIOExceptions = true)
    super.afterAll()
  }

  private def doTheTest(wdlFile: SampleWdl, expectedCalls: TableFor3[String, Long, Long], workflowInputs: Long, workflowOutputs: Long) = {
    val testStart = OffsetDateTime.now
    within(TimeoutDuration) {
      singleWorkflowActor(
        sampleWdl = wdlFile,
        outputFile = Option(metadataFile))
    }
    eventually {
      val metadataFileContent = metadataFile.contentAsString
      val metadata = metadataFileContent.parseJson.asJsObject.fields
      metadata.get("id") shouldNot be(empty)
      metadata.get("status").toStringValue should be("Succeeded")
      metadata.get("submission").toOffsetDateTime should be >= testStart
      val workflowStart = metadata.get("start").toOffsetDateTime
      workflowStart should be >= metadata.get("submission").toOffsetDateTime
      val workflowEnd = metadata.get("end").toOffsetDateTime
      workflowEnd should be >= metadata.get("start").toOffsetDateTime
      metadata.get("inputs").toFields should have size workflowInputs
      metadata.get("outputs").toFields should have size workflowOutputs
      val calls = metadata.get("calls").toFields
      calls should not be empty

      forAll(expectedCalls) { (callName, numInputs, numOutputs) =>
        val callSeq = calls(callName).asInstanceOf[JsArray].elements
        callSeq should have size 1
        val call = callSeq.head.asJsObject.fields
        val inputs = call.get("inputs").toFields
        inputs should have size numInputs
        call.get("executionStatus").toStringValue should be("Done")
        call.get("backend").toStringValue should be("Local")
        call.get("backendStatus").toStringValue should be("Done")
        call.get("outputs").toFields should have size numOutputs
        val callStart = call.get("start").toOffsetDateTime
        callStart should be >= workflowStart
        val callEnd = call.get("end").toOffsetDateTime
        callEnd should be >= callStart
        callEnd should be <= workflowEnd
        call.get("jobId") shouldNot be(empty)
        call("returnCode").asInstanceOf[JsNumber].value should be(0)
        call.get("stdout") shouldNot be(empty)
        call.get("stderr") shouldNot be(empty)
        call("attempt").asInstanceOf[JsNumber].value should be(1)
      }
    }
  }

  "A SingleWorkflowRunnerActor" should {
    // TODO WOM: needs FQNs
    "successfully run a workflow outputting metadata" in {
      val expectedCalls = Table(
        ("callName", "numInputs", "numOutputs"),
        ("three_step.wc", 1L, 1L),
        ("three_step.ps", 0L, 1L),
        ("three_step.cgrep", 2L, 1L))

      doTheTest(ThreeStep, expectedCalls, 1L, 3L)
    }
    "run a workflow outputting metadata with no remaining input expressions" in {
      val expectedCalls = Table(
        ("callName", "numInputs", "numOutputs"),
        ("wf.echo", 1L, 1L),
        ("wf.echo2", 1L, 1L))
      doTheTest(ExpressionsInInputs, expectedCalls, 2L, 2L)
    }
  }
}

class SingleWorkflowRunnerActorWithMetadataOnFailureSpec extends SingleWorkflowRunnerActorSpec {
  val metadataFile = tempFile()

  override protected def afterAll() = {
    metadataFile.delete(swallowIOExceptions = true)
    super.afterAll()
  }

  "A SingleWorkflowRunnerActor" should {
    "fail to run a workflow and still output metadata" in {
      val testStart = OffsetDateTime.now
      within(TimeoutDuration) {
        singleWorkflowActor(sampleWdl = GoodbyeWorld, outputFile = Option(metadataFile))
      }

      val metadata = metadataFile.contentAsString.parseJson.asJsObject.fields
      metadata.get("id") shouldNot be(empty)
      metadata.get("status").toStringValue should be("Failed")
      val workflowStart = metadata.get("start").toOffsetDateTime
      workflowStart should be >= metadata.get("submission").toOffsetDateTime
      val workflowEnd = metadata.get("end").toOffsetDateTime
      workflowEnd should be >= metadata.get("start").toOffsetDateTime
      metadata.get("submission").toOffsetDateTime should be >= testStart
      metadata.get("inputs").toFields should have size 0
      metadata.get("outputs").toFields should have size 0
      val calls = metadata.get("calls").toFields
      calls should not be empty

      val callSeq = calls("wf_goodbye.goodbye").asInstanceOf[JsArray].elements
      callSeq should have size 1
      val call = callSeq.head.asJsObject.fields
      val inputs = call.get("inputs").toFields
      inputs should have size 0
      call.get("executionStatus").toStringValue should be("Failed")
      call.get("backend").toStringValue should be("Local")
      call.get("backendStatus").toStringValue should be("Done")
      call.get("outputs") shouldBe empty
      val callStart = call.get("start").toOffsetDateTime
      callStart should be >= workflowStart
      val callEnd = call.get("end").toOffsetDateTime
      callEnd should be >= callStart
      callEnd should be <= workflowEnd
      call.get("jobId") shouldNot be(empty)
      call("returnCode").asInstanceOf[JsNumber].value shouldNot be (0)
      call.get("stdout") shouldNot be(empty)
      call.get("stderr") shouldNot be(empty)
      call("attempt").asInstanceOf[JsNumber].value should be (1)
      call("failures").asInstanceOf[JsArray].elements shouldNot be(empty)
    }
  }
}

class SingleWorkflowRunnerActorWithBadMetadataSpec extends SingleWorkflowRunnerActorSpec {
  val metadataDir = tempDir()

  override protected def afterAll() = {
    metadataDir.delete(swallowIOExceptions = true)
    super.afterAll()
  }

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow requesting a bad metadata path" in {
      within(TimeoutDuration) {
        val runner = createRunnerActor(outputFile = Option(metadataDir))
        waitForErrorWithException(s"Specified metadata path is a directory, should be a file: $metadataDir") {
          val futureResult = runner.ask(RunWorkflow)(30.seconds.dilated, implicitly)
          Await.ready(futureResult, Duration.Inf)
          futureResult.value.get match {
            case Success(_) =>
            case Failure(e) =>
              e.printStackTrace()
              fail(e)
          }
        }
      }
    }
  }
}

class SingleWorkflowRunnerActorFailureSpec extends SingleWorkflowRunnerActorSpec {
  "A SingleWorkflowRunnerActor" should {
    "successfully terminate the system on an exception" in {
      within(TimeoutDuration) {
        val runner = createRunnerActor()
        val futureResult = runner ? RunWorkflow
        val ex = new RuntimeException("expected error") with NoStackTrace
        runner ! Status.Failure(ex)
        Await.ready(futureResult, Duration.Inf)
        futureResult.value.get match {
          case Success(_) => fail("Unexpected success")
          case Failure(e) => e.getMessage should include("expected error")
        }
      }
    }
  }
}

class SingleWorkflowRunnerActorUnexpectedSpec extends SingleWorkflowRunnerActorSpec {
  "A SingleWorkflowRunnerActor" should {
    "successfully warn about unexpected output" in {
      within(TimeoutDuration) {
        val runner = createRunnerActor()
        waitForWarning("SingleWorkflowRunnerActor: received unexpected message: expected unexpected") {
          runner ? RunWorkflow
          runner ! "expected unexpected"
        }
        assert(!system.whenTerminated.isCompleted)
      }
    }
  }
}
