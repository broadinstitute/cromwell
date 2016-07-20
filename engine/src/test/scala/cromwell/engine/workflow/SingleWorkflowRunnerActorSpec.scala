package cromwell.engine.workflow

import java.nio.file.Path
import java.time.OffsetDateTime

import akka.actor._
import akka.pattern.ask
import akka.testkit.TestKit
import better.files._
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.CromwellTestkitSpec._
import cromwell.engine.workflow.SingleWorkflowRunnerActor.RunWorkflow
import cromwell.engine.workflow.SingleWorkflowRunnerActorSpec._
import cromwell.engine.workflow.workflowstore.{InMemoryWorkflowStore, WorkflowStoreActor}
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.{ExpressionsInInputs, GoodbyeWorld, ThreeStep}
import org.scalatest.prop.{TableDrivenPropertyChecks, TableFor3}
import spray.json._

import scala.concurrent.Await
import scala.concurrent.duration.Duration
import scala.language.postfixOps
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

  def tempFile() = File.newTemp("metadata.", ".json").path

  def tempDir() = File.newTempDir("metadata.dir.").path

  implicit class OptionJsValueEnhancer(val jsValue: Option[JsValue]) extends AnyVal {
    def toOffsetDateTime = OffsetDateTime.parse(jsValue.toStringValue)
    def toStringValue = jsValue.get.asInstanceOf[JsString].value
    def toFields = jsValue.get.asJsObject.fields
  }
}

abstract class SingleWorkflowRunnerActorSpec extends CromwellTestkitSpec {
  private val workflowStore = system.actorOf(WorkflowStoreActor.props(new InMemoryWorkflowStore))
  def workflowManagerActor(): ActorRef = {
    system.actorOf(Props(new WorkflowManagerActor(ConfigFactory.load(), workflowStore)), "WorkflowManagerActor")
  }
  
  def createRunnerActor(sampleWdl: SampleWdl = ThreeStep, managerActor: => ActorRef = workflowManagerActor(),
                          outputFile: => Option[Path] = None): ActorRef = {
    system.actorOf(SingleWorkflowRunnerActor.props(sampleWdl.asWorkflowSources(), outputFile, managerActor, workflowStore))
  }

  def singleWorkflowActor(sampleWdl: SampleWdl = ThreeStep, managerActor: => ActorRef = workflowManagerActor(),
                          outputFile: => Option[Path] = None): Unit = {
    val actorRef = createRunnerActor(sampleWdl, managerActor, outputFile)
    val futureResult = actorRef ? RunWorkflow
    Await.ready(futureResult, Duration.Inf)
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
    metadataFile.delete(ignoreIOExceptions = true)
    super.afterAll()
  }

  private def doTheTest(wdlFile: SampleWdl, expectedCalls: TableFor3[String, Int, Int], workflowInputs: Int, workflowOutputs: Int) = {
    val testStart = OffsetDateTime.now
    within(TimeoutDuration) {
      singleWorkflowActor(
        sampleWdl = wdlFile,
        outputFile = Option(metadataFile))
    }
    TestKit.shutdownActorSystem(system, TimeoutDuration)

    val metadata = metadataFile.contentAsString.parseJson.asJsObject.fields
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
      val callSeq = calls.get(callName).get.asInstanceOf[JsArray].elements
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
      call.get("jobId") should be(empty)
      call.get("returnCode").get.asInstanceOf[JsNumber].value should be (0)
      call.get("stdout") shouldNot be(empty)
      call.get("stderr") shouldNot be(empty)
      call.get("attempt").get.asInstanceOf[JsNumber].value should be (1)
    }
  }

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow outputting metadata" in {
      val expectedCalls = Table(
        ("callName", "numInputs", "numOutputs"),
        ("three_step.wc", 1, 1),
        ("three_step.ps", 0, 1),
        ("three_step.cgrep", 2, 1))

      doTheTest(ThreeStep, expectedCalls, 1, 3)
    }
    "run a workflow outputting metadata with no remaining input expressions" in {
      val expectedCalls = Table(
        ("callName", "numInputs", "numOutputs"),
        ("wf.echo", 1, 1),
        ("wf.echo2", 1, 1))
      doTheTest(ExpressionsInInputs, expectedCalls, 2, 2)
    }
  }
}

class SingleWorkflowRunnerActorWithMetadataOnFailureSpec extends SingleWorkflowRunnerActorSpec {
  val metadataFile = tempFile()

  override protected def afterAll() = {
    metadataFile.delete(ignoreIOExceptions = true)
    super.afterAll()
  }

  "A SingleWorkflowRunnerActor" should {
    "fail to run a workflow and still output metadata" in {
      val testStart = OffsetDateTime.now
      within(TimeoutDuration) {
        singleWorkflowActor(sampleWdl = GoodbyeWorld, outputFile = Option(metadataFile))
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

      val callSeq = calls.get("goodbye.goodbye").get.asInstanceOf[JsArray].elements
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
      call.get("jobId") should be(empty)
      call.get("returnCode").get.asInstanceOf[JsNumber].value shouldNot be (0)
      call.get("stdout") shouldNot be(empty)
      call.get("stderr") shouldNot be(empty)
      call.get("attempt").get.asInstanceOf[JsNumber].value should be (1)
      call.get("failures").get.asInstanceOf[JsArray].elements shouldNot be(empty)
    }
  }
}

class SingleWorkflowRunnerActorWithBadMetadataSpec extends SingleWorkflowRunnerActorSpec {
  val metadataDir = tempDir()

  override protected def afterAll() = {
    metadataDir.delete(ignoreIOExceptions = true)
    super.afterAll()
  }

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow requesting a bad metadata path" in {
      within(TimeoutDuration) {
        val runner = createRunnerActor(outputFile = Option(metadataDir))
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
        assert(!system.isTerminated)
      }
      TestKit.shutdownActorSystem(system, TimeoutDuration)
    }
  }
}
