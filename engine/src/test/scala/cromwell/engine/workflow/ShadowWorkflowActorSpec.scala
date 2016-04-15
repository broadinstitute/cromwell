package cromwell.engine.workflow

import akka.testkit.TestFSMRef
import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.WorkflowSourceFiles
import cromwell.engine.backend.WorkflowDescriptorBuilder
import cromwell.engine.workflow.MaterializeWorkflowDescriptorActor.MaterializeWorkflowDescriptorSuccess
import cromwell.engine.workflow.ShadowWorkflowActor._
import cromwell.util.SampleWdl.{ThreeStep, HelloWorld}

import spray.json._
import spray.json.DefaultJsonProtocol._
import wdl4s.Call

import scala.collection._

class ShadowWorkflowActorSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder {

  override implicit val actorSystem = system

  private def createWorkflowActor(workflowId: WorkflowId = WorkflowId.randomId(),
                                  startMode: StartMode = StartNewWorkflow,
                                  wdlSource: WorkflowSourceFiles) = {
    TestFSMRef(new ShadowWorkflowActor(workflowId, startMode, wdlSource))
  }

  "ShadowWorkflowActor" should {
    "move from WorkflowUnstartedState to MaterializingWorkflowDescriptorState" in {
      val wdlSources = WorkflowSourceFiles(HelloWorld.wdlSource("runtime { }"),
        HelloWorld.rawInputs.toJson.toString(), "{ }")
      val testActor = createWorkflowActor(wdlSource = wdlSources)
      testActor.setState(stateName = WorkflowUnstartedState, stateData = ShadowWorkflowActorData(Map.empty))
      testActor ! StartWorkflowCommand
      testActor.stateName should be (ShadowWorkflowActor.MaterializingWorkflowDescriptorState)
      testActor.stop()
    }
    "be happy to produce correct call assignments if workflow options are defined" in {
      import scala.concurrent.duration._
      val wfOptions =
        """{
          | "backend": "Jes"
          |}
        """.stripMargin
      val runtimeAttributes =
        """runtime {
          | backend: "SGE"
          |}""".stripMargin
      // WorkflowOptions (Jes) should trump runtime attributes (SGE)
      val wdlSources = ThreeStep.asWorkflowSources(runtime = runtimeAttributes, workflowOptions = wfOptions)
      val descriptor = materializeWorkflowDescriptorFromSources(workflowSources = wdlSources)
      val expectedStateData = ShadowWorkflowActorData(descriptor.namespace.workflow.calls.map(_ -> "Jes") (breakOut): scala.collection .immutable.Map[Call, String])
      val testActor = createWorkflowActor(wdlSource = wdlSources)
      testActor.setState(stateName = MaterializingWorkflowDescriptorState, stateData = ShadowWorkflowActorData(Map.empty))
      within(5.seconds) {
        testActor ! MaterializeWorkflowDescriptorSuccess(descriptor)
        testActor.stateName should be(ShadowWorkflowActor.InitializingWorkflowState)
        testActor.stateData should be(expectedStateData)
      }
      testActor.stop()
    }
  }
}
