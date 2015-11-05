package cromwell.engine.workflow

import java.nio.file.Path

import akka.actor._
import akka.testkit._
import better.files._
import cromwell.CromwellTestkitSpec
import cromwell.CromwellTestkitSpec._
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.workflow.SingleWorkflowRunnerActorSpec._
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.{GoodbyeWorld, ThreeStep}
import cromwell.webservice.WorkflowJsonSupport._
import cromwell.webservice.WorkflowMetadataResponse
import org.scalatest.prop.TableDrivenPropertyChecks
import spray.json._

import scala.language.postfixOps

/**
 * A series of tests of the SingleWorkflowRunnerActor. Currently uses live versions of the SingleWorkflowRunnerActor and
 * the WorkflowManagerActor communicating with each other, instead of TestActor/TestProbe.
 *
 * Currently, as instance of the actor system are created via an instance of CromwellTestkitSpec, and the
 * SingleWorkflowRunnerActor also tests halting its actor system, each spec is currently in a separate instance of the
 * CromwellTestkitSpec.
 */
object SingleWorkflowRunnerActorSpec {

  class SilentReceivingActor extends Actor {
    override def receive = {
      case WorkflowManagerActor.Shutdown => context.system.shutdown()
      case _ => /* consume everything else and don't reply */
    }
  }

  def tempFile() = File.newTemp("metadata.", ".json").path

  def tempDir() = File.newTempDir("metadata.dir.").path
}

abstract class SingleWorkflowRunnerActorSpec(name: String) extends CromwellTestkitSpec(name) {
  def workflowManagerActor(): ActorRef = {
    system.actorOf(Props(classOf[WorkflowManagerActor], new LocalBackend))
  }

  def silentActor(): ActorRef = {
    system.actorOf(Props(classOf[SilentReceivingActor]))
  }

  def singleWorkflowActor(sampleWdl: SampleWdl = ThreeStep, managerActor: => ActorRef = workflowManagerActor(),
                          outputFile: => Option[Path] = None): ActorRef = {
    system.actorOf(SingleWorkflowRunnerActor.props(sampleWdl.asWorkflowSources(), outputFile, managerActor))
  }
}

class SingleWorkflowRunnerActorNormalSpec extends SingleWorkflowRunnerActorSpec("SingleWorkflowRunnerActorNormalSpec") {
  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow" in {
      within(timeoutDuration) {
        waitForInfo("workflow finished with status 'Succeeded'.") {
          singleWorkflowActor()
          system.awaitTermination()
        }
      }
    }
  }
}

class SingleWorkflowRunnerActorWithMetadataSpec extends
SingleWorkflowRunnerActorSpec("SingleWorkflowRunnerActorWithMetadataSpec") with TableDrivenPropertyChecks {
  val metadataFile = tempFile()

  override protected def afterAll() = metadataFile.delete(ignoreIOExceptions = true)

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow outputting metadata" in {
      val testStart = System.currentTimeMillis
      within(timeoutDuration) {
        singleWorkflowActor(outputFile = Option(metadataFile))
        system.awaitTermination()
      }

      val metadata = metadataFile.contentAsString.parseJson.convertTo[WorkflowMetadataResponse]
      metadata.id shouldNot be(empty)
      metadata.status should be("Succeeded")
      metadata.submission.getMillis should be >= testStart
      metadata.start shouldNot be(empty)
      metadata.start.get.getMillis should be >= metadata.submission.getMillis
      metadata.end shouldNot be(empty)
      metadata.end.get.getMillis should be >= metadata.start.get.getMillis
      metadata.inputs.fields should have size 1
      metadata.outputs shouldNot be(empty)
      metadata.outputs.get should have size 3
      metadata.calls shouldNot be(empty)

      val expectedCalls = Table(
        ("callName", "numInputs", "numOutputs"),
        ("three_step.wc", 1, 1),
        ("three_step.ps", 0, 1),
        ("three_step.cgrep", 2, 1))

      forAll(expectedCalls) { (callName, numInputs, numOutputs) =>
        val callSeq = metadata.calls(callName)
        callSeq should have size 1
        val call = callSeq.head
        call.inputs should have size numInputs
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
      }
    }
  }
}

class SingleWorkflowRunnerActorWithMetadataOnFailureSpec extends
SingleWorkflowRunnerActorSpec("SingleWorkflowRunnerActorWithMetadataOnFailureSpec") {
  val metadataFile = tempFile()

  override protected def afterAll() = metadataFile.delete(ignoreIOExceptions = true)

  "A SingleWorkflowRunnerActor" should {
    "fail to run a workflow and still output metadata" in {
      val testStart = System.currentTimeMillis
      within(timeoutDuration) {
        singleWorkflowActor(sampleWdl = GoodbyeWorld, outputFile = Option(metadataFile))
        system.awaitTermination()
      }

      val metadata = metadataFile.contentAsString.parseJson.convertTo[WorkflowMetadataResponse]
      metadata.id shouldNot be(empty)
      metadata.status should be("Failed")
      metadata.submission.getMillis should be >= testStart
      metadata.start shouldNot be(empty)
      metadata.start.get.getMillis should be >= metadata.submission.getMillis
      if (metadata.end.isDefined) {
        metadata.end.get.getMillis should be >= metadata.start.get.getMillis
      }
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
      if (metadata.end.isDefined) {
        call.end.get.getMillis should be <= metadata.end.get.getMillis
      }
      call.jobId should be(empty)
      call.returnCode shouldNot be(empty)
      call.returnCode.get shouldNot be(0)
      call.stdout shouldNot be(empty)
      call.stderr shouldNot be(empty)
    }
  }
}

class SingleWorkflowRunnerActorWithBadMetadataSpec extends
SingleWorkflowRunnerActorSpec("SingleWorkflowRunnerActorWithBadMetadataSpec") {
  val metadataDir = tempDir()

  override protected def afterAll() = metadataDir.delete(ignoreIOExceptions = true)

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow requesting a bad metadata path" in {
      within(timeoutDuration) {
        waitForErrorWithException(s"SingleWorkflowRunnerActor: $metadataDir: Is a directory") {
          singleWorkflowActor(outputFile = Option(metadataDir))
          system.awaitTermination()
        }
      }
    }
  }
}

class SingleWorkflowRunnerActorFailureSpec extends
SingleWorkflowRunnerActorSpec("SingleWorkflowRunnerActorFailureSpec") {
  "A SingleWorkflowRunnerActor" should {
    "successfully terminate the system on an exception" in {
      within(timeoutDuration) {
        val actor = singleWorkflowActor(managerActor = silentActor())
        val ex = new RuntimeException("expected error")
        ex.setStackTrace(ex.getStackTrace.take(1)) // Leave just a small hint, not a full trace.
        actor ! Status.Failure(ex)
        system.awaitTermination()
      }
    }
  }
}

class SingleWorkflowRunnerActorUnexpectedSpec extends
SingleWorkflowRunnerActorSpec("SingleWorkflowRunnerActorUnexpectedSpec") {
  "A SingleWorkflowRunnerActor" should {
    "successfully warn about unexpected output" in {
      within(timeoutDuration) {
        waitForWarning("SingleWorkflowRunnerActor: received unexpected message: expected unexpected") {
          val actor = singleWorkflowActor(managerActor = silentActor())
          actor ! "expected unexpected"
        }
        assert(!system.isTerminated)
      }
      TestKit.shutdownActorSystem(system, timeoutDuration)
    }
  }
}
