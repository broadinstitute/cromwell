package cromwell.engine.workflow

import java.nio.file.Path
import java.time.OffsetDateTime

import akka.actor._
import akka.pattern.ask
import akka.testkit.TestKit
import akka.util.Timeout
import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestKitSpec._
import cromwell.core.{WorkflowSourceFilesCollection}
import cromwell.engine.backend.BackendSingletonCollection
import cromwell.engine.workflow.SingleWorkflowRunnerActor.RunWorkflow
import cromwell.engine.workflow.SingleWorkflowRunnerActorSpec._
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor
import cromwell.engine.workflow.workflowstore.{InMemoryWorkflowStore, WorkflowStoreActor}
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.{ExpressionsInInputs, GoodbyeWorld, ThreeStep}
import cromwell.{AlwaysHappyJobStoreActor, AlwaysHappySubWorkflowStoreActor, CromwellTestKitSpec, EmptyCallCacheReadActor}
import org.scalatest.prop.{TableDrivenPropertyChecks, TableFor3}
import spray.json._

import scala.concurrent.Await
import scala.concurrent.duration.Duration
import scala.util._

/**
 * A series of tests of the SingleWorkflowRunnerActor. Currently uses live versions of the SingleWorkflowRunnerActor and
 * the WorkflowManagerActor communicating with each other, instead of TestActor/TestProbe.
 *
 * Currently, as instance of the actor system are created via an instance of CromwellTestkitSpec, and the
 * SingleWorkflowRunnerActor also tests halting its actor system, each spec is currently in a separate instance of the
 * CromwellTestkitSpec.
 */
object SingleWorkflowRunnerActorSpec {

  def tempFile() = File.newTemporaryFile("metadata.", ".json")

  def tempDir() = File.newTemporaryDirectory("metadata.dir.")

  implicit class OptionJsValueEnhancer(val jsValue: Option[JsValue]) extends AnyVal {
    def toOffsetDateTime = OffsetDateTime.parse(jsValue.toStringValue)
    def toStringValue = jsValue.get.asInstanceOf[JsString].value
    def toFields = jsValue.get.asJsObject.fields
  }

  class TestSingleWorkflowRunnerActor(source: WorkflowSourceFilesCollection,
                                      metadataOutputPath: Option[Path])
    extends SingleWorkflowRunnerActor(source, metadataOutputPath) {
    override lazy val serviceRegistryActor = CromwellTestKitSpec.ServiceRegistryActorInstance
  }
}

abstract class SingleWorkflowRunnerActorSpec extends CromwellTestKitSpec {
  private val workflowStore = system.actorOf(WorkflowStoreActor.props(new InMemoryWorkflowStore, dummyServiceRegistryActor))
  private val jobStore = system.actorOf(AlwaysHappyJobStoreActor.props)
  private val subWorkflowStore = system.actorOf(AlwaysHappySubWorkflowStoreActor.props)
  private val callCacheReadActor = system.actorOf(EmptyCallCacheReadActor.props)
  private val jobTokenDispenserActor = system.actorOf(JobExecutionTokenDispenserActor.props)


  def workflowManagerActor(): ActorRef = {
    val params = WorkflowManagerActorParams(ConfigFactory.load(),
      workflowStore,
      dummyServiceRegistryActor,
      dummyLogCopyRouter,
      jobStore,
      subWorkflowStore,
      callCacheReadActor,
      jobTokenDispenserActor,
      BackendSingletonCollection(Map.empty),
      abortJobsOnTerminate = false,
      serverMode = false)
    system.actorOf(Props(new WorkflowManagerActor(params)), "WorkflowManagerActor")
  }
  
  def createRunnerActor(sampleWdl: SampleWdl = ThreeStep, managerActor: => ActorRef = workflowManagerActor(),
                          outputFile: => Option[Path] = None): ActorRef = {
    system.actorOf(Props(new TestSingleWorkflowRunnerActor(sampleWdl.asWorkflowSources(), outputFile)))
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
      TestKit.shutdownActorSystem(system, TimeoutDuration)
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
        outputFile = Option(metadataFile.path))
        TestKit.shutdownActorSystem(system, TimeoutDuration)
    }
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
      call.get("backendStatus") should be(empty)
      call.get("outputs").toFields should have size numOutputs
      val callStart = call.get("start").toOffsetDateTime
      callStart should be >= workflowStart
      val callEnd = call.get("end").toOffsetDateTime
      callEnd should be >= callStart
      callEnd should be <= workflowEnd
      call.get("jobId") shouldNot be(empty)
      call("returnCode").asInstanceOf[JsNumber].value should be (0)
      call.get("stdout") shouldNot be(empty)
      call.get("stderr") shouldNot be(empty)
      call("attempt").asInstanceOf[JsNumber].value should be (1)
    }
  }

  "A SingleWorkflowRunnerActor" should {
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
        singleWorkflowActor(sampleWdl = GoodbyeWorld, outputFile = Option(metadataFile.path))
      }
      TestKit.shutdownActorSystem(system, TimeoutDuration)

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
      call.get("backendStatus") should be(empty)
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
        val runner = createRunnerActor(outputFile = Option(metadataDir.path))
        waitForErrorWithException(s"Specified metadata path is a directory, should be a file: $metadataDir") {
          val futureResult = runner ? RunWorkflow
          Await.ready(futureResult, Duration.Inf)
          futureResult.value.get match {
            case Success(_) =>
            case Failure(e) =>
              e.printStackTrace()
              fail(e)
          }
        }
      }
      TestKit.shutdownActorSystem(system, TimeoutDuration)
    }
  }
}

class SingleWorkflowRunnerActorFailureSpec extends SingleWorkflowRunnerActorSpec {
  "A SingleWorkflowRunnerActor" should {
    "successfully terminate the system on an exception" in {
      within(TimeoutDuration) {
        val runner = createRunnerActor()
        val futureResult = runner ? RunWorkflow
        val ex = new RuntimeException("expected error")
        ex.setStackTrace(ex.getStackTrace.take(1)) // Leave just a small hint, not a full trace.
        runner ! Status.Failure(ex)
        Await.ready(futureResult, Duration.Inf)
        futureResult.value.get match {
          case Success(_) => fail("Unexpected success")
          case Failure(e) => e.getMessage should include("expected error")
        }
      }
      TestKit.shutdownActorSystem(system, TimeoutDuration)
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
      TestKit.shutdownActorSystem(system, TimeoutDuration)
    }
  }
}
