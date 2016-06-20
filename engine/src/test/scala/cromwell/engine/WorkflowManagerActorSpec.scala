package cromwell.engine

import java.util.UUID

import akka.testkit._
import cromwell.CromwellSpec.DockerTest
import cromwell.CromwellTestkitSpec
import cromwell.CromwellTestkitSpec._
import cromwell.core.ExecutionStatus.{NotStarted, Running}
import cromwell.core._
import cromwell.engine.backend.WorkflowDescriptorBuilder
import cromwell.engine.db.DataAccess._
import cromwell.util.SampleWdl
import cromwell.util.SampleWdl.{HelloWorld, HelloWorldWithoutWorkflow, Incr}
import cromwell.webservice.CromwellApiHandler._
import wdl4s._
import wdl4s.types.{WdlArrayType, WdlStringType}
import wdl4s.values.{SymbolHash, WdlArray, WdlInteger, WdlString}

import scala.concurrent.duration._
import scala.concurrent.{Await, Future}
import scala.language.postfixOps

class WorkflowManagerActorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder {
  override implicit val actorSystem = system

  "A WorkflowManagerActor" should {

    val TestExecutionTimeout = 5.seconds.dilated

    "run the Hello World workflow" ignore {
//
//      implicit val workflowManagerActor = TestActorRef(OldStyleWorkflowManagerActor.props(), self, "Test the WorkflowManagerActor")
//
//      val workflowId = waitForHandledMessagePattern(pattern = "transitioning from Running to Succeeded") {
//        messageAndWait[WorkflowManagerSubmitSuccess](SubmitWorkflow(HelloWorld.asWorkflowSources())).id
//      }
//
//      val status = messageAndWait[WorkflowManagerStatusSuccess](WorkflowStatus(workflowId)).state
//      status shouldEqual WorkflowSucceeded
//
//      val workflowOutputs = messageAndWait[WorkflowManagerWorkflowOutputsSuccess](WorkflowOutputs(workflowId)).outputs
//
//      val actualWorkflowOutputs = workflowOutputs.map { case (k, JobOutput(WdlString(string), _)) => k -> string }
//      actualWorkflowOutputs shouldEqual Map(HelloWorld.OutputKey -> HelloWorld.OutputValue)
//
//      val callOutputs = messageAndWait[WorkflowManagerCallOutputsSuccess](CallOutputs(workflowId, "hello.hello")).outputs
//      val actualCallOutputs = callOutputs.map { case (k, JobOutput(WdlString(string), _)) => k -> string }
//      actualCallOutputs shouldEqual Map("salutation" -> HelloWorld.OutputValue)
    }

    "Not try to restart any workflows when there are no workflows in restartable states" ignore {
//      waitForPattern("Found no workflows to restart.") {
//        TestActorRef(OldStyleWorkflowManagerActor.props(), self, "No workflows")
//      }
    }

    "Try to restart workflows when there are workflows in restartable states" ignore {
//      val workflows = Map(
//        WorkflowId(UUID.randomUUID()) -> WorkflowSubmitted,
//        WorkflowId(UUID.randomUUID()) -> WorkflowRunning)
//      val ids = workflows.keys.map(_.toString).toSeq.sorted
//      val key = SymbolStoreKey("hello.hello", "addressee", None, input = true)
//      val worldWdlString = WdlString("world")
//      val setupFuture = Future.sequence(
//        workflows map { case (workflowId, workflowState) =>
//          val status = if (workflowState == WorkflowSubmitted) NotStarted else Running
//          val sources = SampleWdl.HelloWorld.asWorkflowSources()
//          val descriptor = createMaterializedEngineWorkflowDescriptor(id = workflowId, workflowSources = sources)
//          // PBE hacked out
//          // val worldSymbolHash = descriptor.hash(worldWdlString)
//          val worldSymbolHash = Option(SymbolHash(UUID.randomUUID().toString))
//          val symbols = Map(key -> new SymbolStoreEntry(key, WdlStringType, Option(worldWdlString), worldSymbolHash))
//          // FIXME? null AST
//          val task = Task.empty
//          val call = new Call(None, key.scope, task, Set.empty[FullyQualifiedName], Map.empty, None)
//          for {
//            _ <- globalDataAccess.createWorkflow(descriptor, sources, symbols.values, Seq(call), new OldStyleLocalBackend(CromwellTestkitSpec.DefaultLocalBackendConfigEntry, system))
//            _ <- globalDataAccess.updateWorkflowState(workflowId, workflowState)
//            _ <- globalDataAccess.updateStatus(workflowId, Seq(ExecutionDatabaseKey(call.fullyQualifiedName, None, 1)), status)
//          } yield ()
//        }
//      )
//      Await.result(setupFuture, Duration.Inf)
//
//      waitForPattern("Restarting workflow IDs: " + ids.mkString(", ")) {
//        waitForPattern("Found 2 workflows to restart.") {
//          // Workflows are always set back to Submitted on restart.
//          waitForPattern("transitioning from Submitted to Running.", occurrences = 2) {
//            // Both the previously in-flight call and the never-started call should get started.
//            waitForPattern("starting calls: hello.hello", occurrences = 2) {
//              waitForPattern("transitioning from Running to Succeeded", occurrences = 2) {
//                TestActorRef(OldStyleWorkflowManagerActor.props(), self, "2 restartable workflows")
//              }
//            }
//          }
//        }
//      }
    }


    "Handle coercion failures gracefully" ignore {
//      within(TestExecutionTimeout) {
//        implicit val workflowManagerActor = TestActorRef(OldStyleWorkflowManagerActor.props(), self, "Test WorkflowManagerActor coercion failures")
//        waitForErrorWithException("Workflow failed submission") {
//          val e = messageAndWait[WorkflowManagerSubmitFailure](SubmitWorkflow(Incr.asWorkflowSources())).failure
//          e.getMessage should include("Could not coerce value for 'incr.incr.val' into: WdlIntegerType")
//        }
//      }
    }

    "error when running a workflowless WDL" ignore {
//      implicit val workflowManagerActor = TestActorRef(OldStyleWorkflowManagerActor.props(), self, "Test a workflowless submission")
//      val e = messageAndWait[WorkflowManagerSubmitFailure](SubmitWorkflow(HelloWorldWithoutWorkflow.asWorkflowSources())).failure
//      e.getMessage should include("Namespace does not have a local workflow to run")
    }

    "error when asked for outputs of a nonexistent workflow" ignore {
//      within(TestExecutionTimeout) {
//        implicit val workflowManagerActor = TestActorRef(OldStyleWorkflowManagerActor.props(),
//          self, "Test WorkflowManagerActor output lookup failure")
//        val id = WorkflowId(UUID.randomUUID())
//
//        val e = messageAndWait[WorkflowManagerWorkflowOutputsFailure](WorkflowOutputs(id)).failure
//        e.getMessage shouldBe s"Workflow '$id' not found"
//      }
    }

    "error when asked for call logs of a nonexistent workflow" ignore {
//      within(TestExecutionTimeout) {
//        implicit val workflowManagerActor = TestActorRef(OldStyleWorkflowManagerActor.props(),
//          self, "Test WorkflowManagerActor call log lookup failure")
//        val id = WorkflowId.randomId()
//        val e = messageAndWait[WorkflowManagerCallStdoutStderrFailure](CallStdoutStderr(id, "foo.bar")).failure
//        e.getMessage shouldBe s"Workflow '$id' not found"
//      }
    }


    "error when asked for logs of a nonexistent workflow" ignore {
//      within(TestExecutionTimeout) {
//        implicit val workflowManagerActor = TestActorRef(OldStyleWorkflowManagerActor.props(),
//          self, "Test WorkflowManagerActor log lookup failure")
//        val id = WorkflowId.randomId()
//        val e = messageAndWait[WorkflowManagerWorkflowStdoutStderrFailure](WorkflowStdoutStderr(id)).failure
//        e.getMessage shouldBe s"Workflow '$id' not found"
//      }
    }

    "run workflows in the correct directory" ignore {
//      val outputs = runWdl(sampleWdl = SampleWdl.CurrentDirectory,
//        EventFilter.info(pattern = s"starting calls: whereami.whereami", occurrences = 1))
//      val outputName = "whereami.whereami.pwd"
//      val salutation = outputs.get(outputName).get
//      val actualOutput = salutation.asInstanceOf[JobOutput].wdlValue.asInstanceOf[WdlString].value.trim
//      actualOutput should endWith("/call-whereami")
    }

    "build metadata correctly" ignore {
//
//      implicit val workflowManagerActor = TestActorRef(OldStyleWorkflowManagerActor.props(), self, "Test Workflow metadata construction")
//
//      val workflowId = waitForHandledMessagePattern(pattern = "transitioning from Running to Succeeded") {
//        messageAndWait[WorkflowManagerSubmitSuccess](SubmitWorkflow(new SampleWdl.ScatterWdl().asWorkflowSources())).id
//      }
//
//      val status = messageAndWait[WorkflowManagerStatusSuccess](WorkflowStatus(workflowId)).state
//      status shouldEqual WorkflowSucceeded
//
//      val metadata = messageAndWait[WorkflowManagerWorkflowMetadataSuccess](WorkflowMetadata(workflowId)).response
//      metadata should not be null
//
//      metadata.status shouldBe WorkflowSucceeded.toString
//      metadata.start shouldBe defined
//      metadata.end shouldBe defined
//      metadata.outputs shouldBe defined
//      metadata.outputs.get should have size 5
//      metadata.calls should have size 5
//
//      // ok if this explodes, it's a test
//      val devOutputs = metadata.outputs.get.get("w.A.A_out").get
//      val wdlArray = devOutputs.asInstanceOf[WdlArray]
//      wdlArray.wdlType shouldBe WdlArrayType(WdlStringType)
//
//      (wdlArray.value map { case WdlString(string) => string }) shouldEqual Vector("jeff", "chris", "miguel", "thibault", "khalid", "scott")
//
//      val devCalls = metadata.calls.get("w.C").get
//      devCalls should have size 6
//      devCalls foreach { call =>
//        call.start shouldBe defined
//        call.end shouldBe defined
//        call.jobId should not be defined
//        call.returnCode.get shouldBe 0
//        call.stdout shouldBe defined
//        call.stderr shouldBe defined
//        call.inputs should have size 1
//        call.backend.get shouldEqual "Local"
//        call.backendStatus should not be defined
//        call.executionStatus shouldBe "Done"
//        call.attempt shouldBe 1
//      }
//
//      (devCalls map { _.outputs.get.get("C_out").get.asInstanceOf[WdlInteger].value }) shouldEqual Vector(400, 500, 600, 800, 600, 500)
    }

    "show (only supported) runtime attributes in metadata" taggedAs DockerTest ignore {
//
//      implicit val workflowManagerActor = TestActorRef(OldStyleWorkflowManagerActor.props(), self, "Test Workflow metadata construction")
//
//      val fullWfOptions =
//        """
//          |{
//          |  "default_runtime_attributes": {
//          |    "continueOnReturnCode": [0, 1, 2, 3],
//          |    "failOnStderr": true,
//          |    "zones": "us-central1-a",
//          |    "disks": "local-disk 10 SSD",
//          |    "memory": "5000000 KB",
//          |    "preemptible": 2,
//          |    "cpu": 3
//          |  }
//          |}
//        """.stripMargin
//
//      val workflowId = waitForHandledMessagePattern(pattern = "transitioning from Running to Succeeded") {
//        messageAndWait[WorkflowManagerSubmitSuccess](SubmitWorkflow(SampleWdl.WorkflowWithStaticRuntime.asWorkflowSources(runtime = "", workflowOptions = fullWfOptions))).id
//      }
//
//      val status = messageAndWait[WorkflowManagerStatusSuccess](WorkflowStatus(workflowId)).state
//      status shouldEqual WorkflowSucceeded
//
//      val metadata = messageAndWait[WorkflowManagerWorkflowMetadataSuccess](WorkflowMetadata(workflowId)).response
//      metadata should not be null
//
//      metadata.status shouldBe WorkflowSucceeded.toString
//      val cgrep: OldStyleCallMetadata = metadata.calls("two_step.cgrep").head
//      cgrep.runtimeAttributes.size shouldBe 3
//      cgrep.runtimeAttributes("continueOnReturnCode") shouldBe "[0, 1, 2, 3]"
//      cgrep.runtimeAttributes("failOnStderr") shouldBe "true"
//      cgrep.runtimeAttributes("docker") shouldBe "ubuntu:latest"
    }
  }
}
