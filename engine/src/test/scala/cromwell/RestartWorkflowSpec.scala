package cromwell

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.core.Tags._
import cromwell.core._
import cromwell.engine.workflow.WorkflowDescriptorBuilder

class RestartWorkflowSpec extends CromwellTestKitSpec with WorkflowDescriptorBuilder {

  val actorSystem = ActorSystem("RestartWorkflowSpec", ConfigFactory.parseString(CromwellTestKitSpec.ConfigText))
  //val localBackend = new OldStyleLocalBackend(CromwellTestkitSpec.DefaultLocalBackendConfigEntry, actorSystem)
  val sources = WorkflowSourceFilesWithoutImports(
    wdlSource="""task a {command{}}
                |workflow w {
                |  call a
                |  call a as b
                |}
              """.stripMargin,
    inputsJson="{}",
    workflowOptionsJson="{}"
  )

  "RestartWorkflowSpec" should {
    "restart a call in Running state" taggedAs PostMVP ignore {
//      val id = WorkflowId.randomId()
//      val descriptor = createMaterializedEngineWorkflowDescriptor(id, sources)
//      val a = ExecutionDatabaseKey("w.a", Option(-1), 1)
//      val b = ExecutionDatabaseKey("w.b", Option(-1), 1)
//
//      (for {
//        _ <- dataAccess.createWorkflow(descriptor, sources, Nil, descriptor.namespace.workflow.calls, localBackend)
//        _ <- dataAccess.updateWorkflowState(descriptor.id, WorkflowRunning)
//        _ <- dataAccess.updateStatus(descriptor.id, Seq(a), ExecutionStatus.Running)
//        _ <- dataAccess.updateStatus(descriptor.id, Seq(b), ExecutionStatus.NotStarted)
//        wma = (new TestWorkflowManagerSystem).workflowManagerActor
//        _ = verifyWorkflowState(wma, descriptor.id, WorkflowSucceeded)
//        aStatus <- dataAccess.getExecutionStatus(descriptor.id, a)
//        _ = aStatus.map(_.executionStatus) shouldEqual Some(ExecutionStatus.Done)
//        bStatus <- dataAccess.getExecutionStatus(descriptor.id, b)
//        _ = bStatus.map(_.executionStatus) shouldEqual Some(ExecutionStatus.Done)
//      } yield ()).futureValue
    }
  }
}
