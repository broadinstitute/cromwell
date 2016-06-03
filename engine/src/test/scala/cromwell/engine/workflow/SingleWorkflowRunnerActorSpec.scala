package cromwell.engine.workflow

import java.nio.file.Path

import akka.actor._
import akka.pattern.ask
import akka.testkit.{TestActorRef, TestKit}
import better.files._
import cromwell.CromwellTestkitSpec
import cromwell.CromwellTestkitSpec._
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.workflow.SingleWorkflowRunnerActor.RunWorkflow
import cromwell.engine.workflow.SingleWorkflowRunnerActorSpec._
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.{ExpressionsInInputs, GoodbyeWorld, ThreeStep}
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.WorkflowMetadataResponse
import org.scalatest.prop.{TableFor3, TableDrivenPropertyChecks}
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
}

abstract class SingleWorkflowRunnerActorSpec extends CromwellTestkitSpec {

  val isServerMode = false

  def workflowManagerActor(): ActorRef = {
    system.actorOf(WorkflowManagerActor.props(isServerMode))
  }
  
  def createRunnerActor(sampleWdl: SampleWdl = ThreeStep, managerActor: => ActorRef = workflowManagerActor(),
                          outputFile: => Option[Path] = None): ActorRef = {
    system.actorOf(SingleWorkflowRunnerActor.props(sampleWdl.asWorkflowSources(), outputFile, managerActor))
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
    "successfully run a workflow" ignore {
      within(timeoutDuration) {
        waitForInfo("workflow finished with status 'Succeeded'.") {
          implicit val workflowManagerActor = TestActorRef(WorkflowManagerActor.props(isServerMode), self, "Test the SingleWorkflowRunnerActor")
        }
      }
      TestKit.shutdownActorSystem(system, timeoutDuration)
    }
  }
}

class SingleWorkflowRunnerActorWithMetadataSpec extends SingleWorkflowRunnerActorSpec with TableDrivenPropertyChecks {
  val metadataFile = tempFile()

  override protected def afterAll() = metadataFile.delete(ignoreIOExceptions = true)

  private def doTheTest(wdlFile: SampleWdl, expectedCalls: TableFor3[String, Int, Int], workflowInputs: Int, workflowOutputs: Int) = {
    val testStart = System.currentTimeMillis
    within(timeoutDuration) {
      singleWorkflowActor(
        sampleWdl = wdlFile,
        outputFile = Option(metadataFile))
    }
    TestKit.shutdownActorSystem(system, timeoutDuration)

    val metadata = metadataFile.contentAsString.parseJson.convertTo[WorkflowMetadataResponse]
    metadata.id shouldNot be(empty)
    metadata.status should be("Succeeded")
    metadata.submission.getMillis should be >= testStart
    metadata.start shouldNot be(empty)
    metadata.start.get.getMillis should be >= metadata.submission.getMillis
    metadata.end shouldNot be(empty)
    metadata.end.get.getMillis should be >= metadata.start.get.getMillis
    metadata.inputs.fields should have size workflowInputs
    metadata.outputs shouldNot be(empty)
    metadata.outputs.get should have size workflowOutputs
    metadata.calls shouldNot be(empty)

    forAll(expectedCalls) { (callName, numInputs, numOutputs) =>
      val callSeq = metadata.calls(callName)
      callSeq should have size 1
      val call = callSeq.head
      call.inputs should have size numInputs
      call.inputs foreach { case (name, value) =>
        value.wdlType.toWdlString should not be "Expression"
      }
      call.executionStatus should be("Done")
      call.backend should be(Option("Local"))
      call.backendStatus should be(empty)
      call.outputs shouldNot be(empty)
      call.outputs.get should have size numOutputs
      call.start shouldNot be(empty)
      call.start.get.getMillis should be >= metadata.start.get.getMillis
      call.end shouldNot be(empty)
      call.end.get.getMillis should be >= call.start.get.getMillis
      call.end.get.getMillis should be <= metadata.end.get.getMillis
      call.jobId should be(empty)
      call.returnCode should be(Option(0))
      call.stdout shouldNot be(empty)
      call.stderr shouldNot be(empty)
      call.attempt shouldBe 1
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
    "run a workflow outputting metadata with no remaining input expressions" ignore {
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

  override protected def afterAll() = metadataFile.delete(ignoreIOExceptions = true)

  "A SingleWorkflowRunnerActor" should {
    "fail to run a workflow and still output metadata" ignore {
      val testStart = System.currentTimeMillis
      within(timeoutDuration) {
        singleWorkflowActor(sampleWdl = GoodbyeWorld, outputFile = Option(metadataFile))
      }
      TestKit.shutdownActorSystem(system, timeoutDuration)

      val metadata = metadataFile.contentAsString.parseJson.convertTo[WorkflowMetadataResponse]
      metadata.id shouldNot be(empty)
      metadata.status should be("Failed")
      metadata.submission.getMillis should be >= testStart
      metadata.start shouldNot be(empty)
      metadata.start.get.getMillis should be >= metadata.submission.getMillis
      metadata.end.get.getMillis should be >= metadata.start.get.getMillis
      metadata.inputs.fields should have size 0
      metadata.outputs shouldNot be(empty)
      metadata.outputs.get should have size 0
      metadata.calls should have size 1

      val callSeq = metadata.calls("goodbye.goodbye")
      callSeq should have size 1
      val call = callSeq.head
      call.inputs should have size 0
      call.executionStatus should be("Failed")
      call.backend should be(Option("Local"))
      call.backendStatus should be(empty)
      call.outputs shouldNot be(empty)
      call.outputs.get should have size 0
      call.start shouldNot be(empty)
      call.start.get.getMillis should be >= metadata.start.get.getMillis
      call.end shouldNot be(empty)
      call.end.get.getMillis should be >= call.start.get.getMillis
      call.end.get.getMillis should be <= metadata.end.get.getMillis
      call.jobId should be(empty)
      call.returnCode shouldNot be(empty)
      call.returnCode.get shouldNot be(0)
      call.stdout shouldNot be(empty)
      call.stderr shouldNot be(empty)
      call.attempt shouldBe 1
      call.failures shouldNot be(empty)
    }
  }
}

class SingleWorkflowRunnerActorWithBadMetadataSpec extends SingleWorkflowRunnerActorSpec {
  val metadataDir = tempDir()

  override protected def afterAll() = metadataDir.delete(ignoreIOExceptions = true)

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow requesting a bad metadata path" ignore {
      within(timeoutDuration) {
        val runner = createRunnerActor(outputFile = Option(metadataDir))
        waitForErrorWithException(s"$metadataDir: Is a directory") {
          val futureResult = runner ? RunWorkflow
          Await.ready(futureResult, Duration.Inf)
          futureResult.value.get match {
            case Success(_) =>
            case Failure(e) => fail(e)
          }
        }
      }
      TestKit.shutdownActorSystem(system, timeoutDuration)
    }
  }
}

class SingleWorkflowRunnerActorFailureSpec extends SingleWorkflowRunnerActorSpec {
  "A SingleWorkflowRunnerActor" should {
    "successfully terminate the system on an exception" ignore {
      within(timeoutDuration) {
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
      TestKit.shutdownActorSystem(system, timeoutDuration)
    }
  }
}

class SingleWorkflowRunnerActorUnexpectedSpec extends SingleWorkflowRunnerActorSpec {
  "A SingleWorkflowRunnerActor" should {
    "successfully warn about unexpected output" ignore {
      within(timeoutDuration) {
        val runner = createRunnerActor()
        waitForWarning("SingleWorkflowRunnerActor: received unexpected message: expected unexpected") {
          val futureResult = runner ? RunWorkflow
          runner ! "expected unexpected"
        }
        assert(!system.isTerminated)
      }
      TestKit.shutdownActorSystem(system, timeoutDuration)
    }
  }
}