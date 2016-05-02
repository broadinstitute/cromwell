package cromwell.engine.workflow

import akka.testkit.TestFSMRef
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec
import cromwell.core.WorkflowId
import cromwell.engine.workflow.ShadowWorkflowActor._
import cromwell.engine.workflow.lifecycle.ShadowMaterializeWorkflowDescriptorActor
import cromwell.engine.workflow.lifecycle.ShadowMaterializeWorkflowDescriptorActor.{ShadowMaterializeWorkflowDescriptorCommand, ShadowMaterializeWorkflowDescriptorFailureResponse, ShadowMaterializeWorkflowDescriptorSuccessResponse, ShadowWorkflowDescriptorMaterializationResult}
import cromwell.engine.workflow.lifecycle.WorkflowLifecycleActor.WorkflowLifecycleActorData
import cromwell.engine.{EngineWorkflowDescriptor, WorkflowSourceFiles}
import cromwell.util.SampleWdl.{HelloWorld, ThreeStep}
import spray.json.DefaultJsonProtocol._
import spray.json._

import scala.concurrent.Await

class ShadowWorkflowActorSpec extends CromwellTestkitSpec {

  // Blocking call to get the `EngineWorkflowDescriptor`
  private def createMaterializedEngineWorkflowDescriptor(workflowSources: WorkflowSourceFiles): EngineWorkflowDescriptor = {
    import akka.pattern.ask

    import scala.concurrent.duration._

    val actor = system.actorOf(ShadowMaterializeWorkflowDescriptorActor.props())
    val workflowDescriptorFuture =
      (actor ? ShadowMaterializeWorkflowDescriptorCommand(WorkflowId.randomId(), workflowSources, ConfigFactory.load))
        .mapTo[ShadowWorkflowDescriptorMaterializationResult]

    Await.result(workflowDescriptorFuture map {
      case ShadowMaterializeWorkflowDescriptorSuccessResponse(workflowDescriptor) => workflowDescriptor
      case ShadowMaterializeWorkflowDescriptorFailureResponse(reason) => throw reason
    }, 3.seconds)
  }

  private def createWorkflowActor(workflowId: WorkflowId = WorkflowId.randomId(),
                                  startMode: StartMode = StartNewWorkflow,
                                  wdlSource: WorkflowSourceFiles) = {
    TestFSMRef(new ShadowWorkflowActor(workflowId, startMode, wdlSource,ConfigFactory.load()))
  }

  "ShadowWorkflowActor" should {
    "move from WorkflowUnstartedState to MaterializingWorkflowDescriptorState" in {
      val wdlSources = WorkflowSourceFiles(HelloWorld.wdlSource("runtime { }"),
        HelloWorld.rawInputs.toJson.toString(), "{ }")
      val testActor = createWorkflowActor(wdlSource = wdlSources)
      testActor.setState(stateName = WorkflowUnstartedState, stateData = ShadowWorkflowActorData.empty)
      testActor ! StartWorkflowCommand
      testActor.stateName should be (ShadowWorkflowActor.MaterializingWorkflowDescriptorState)
      testActor.stop()
    }

    "transition to InitializingWorkflowState with correct call assignments given workflow options" in {
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
      val descriptor = createMaterializedEngineWorkflowDescriptor(workflowSources = wdlSources)
      val testActor = createWorkflowActor(wdlSource = wdlSources)
      testActor.setState(stateName = MaterializingWorkflowDescriptorState, stateData = ShadowWorkflowActorData.empty)
      within(5.seconds) {
        testActor ! ShadowMaterializeWorkflowDescriptorSuccessResponse(descriptor)
        testActor.stateName should be(ShadowWorkflowActor.InitializingWorkflowState)
        testActor.stateData.workflowDescriptor should be(Some(descriptor))
        testActor.stateData.currentLifecycleStateActor.isDefined should be(true)
      }
      testActor.stop()
    }

    "transition to InitializingWorkflowState with correct call assignments given runtime-attributes" in {
      import scala.concurrent.duration._
      val runtimeAttributes =
        """runtime {
          | backend: "SGE"
          |}""".stripMargin
      val wdlSources = ThreeStep.asWorkflowSources(runtime = runtimeAttributes)
      val descriptor = createMaterializedEngineWorkflowDescriptor(workflowSources = wdlSources)
      val testActor = createWorkflowActor(wdlSource = wdlSources)
      testActor.setState(stateName = MaterializingWorkflowDescriptorState, stateData = ShadowWorkflowActorData.empty)
      within(5.seconds) {
        testActor ! ShadowMaterializeWorkflowDescriptorSuccessResponse(descriptor)
        testActor.stateName should be(ShadowWorkflowActor.InitializingWorkflowState)
        testActor.stateData.workflowDescriptor should be(Some(descriptor))
        testActor.stateData.currentLifecycleStateActor.isDefined should be(true)
      }
      testActor.stop()
    }
  }
}
