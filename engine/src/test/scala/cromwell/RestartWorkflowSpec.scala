package cromwell

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory
import cromwell.CromwellTestkitSpec.TestWorkflowManagerSystem
import cromwell.core.WorkflowId
import cromwell.engine._
import cromwell.engine.backend.WorkflowDescriptorBuilder
import cromwell.engine.backend.local.LocalBackend
import cromwell.engine.db.{DataAccess, ExecutionDatabaseKey}

class RestartWorkflowSpec extends CromwellTestkitSpec with WorkflowDescriptorBuilder {

  val actorSystem = ActorSystem("RestartWorkflowSpec", ConfigFactory.parseString(CromwellTestkitSpec.ConfigText))
  val dataAccess = DataAccess.globalDataAccess
  val localBackend = new LocalBackend(CromwellTestkitSpec.DefaultLocalBackendConfigEntry, actorSystem)
  val sources = WorkflowSourceFiles(
    wdlSource="""task a {command{}}
                |workflow w {
                |  call a
                |  call a as b
                |}
              """.stripMargin,
    inputsJson="{}",
    workflowOptionsJson="{}"
  )

  override protected def afterAll() = {
    actorSystem.shutdown()
    super.afterAll()
  }

  "RestartWorkflowSpec" should {
    "restart a call in Running state" in {
      val id = WorkflowId.randomId()
      val descriptor = materializeWorkflowDescriptorFromSources(id, sources)
      val a = ExecutionDatabaseKey("w.a", Option(-1), 1)
      val b = ExecutionDatabaseKey("w.b", Option(-1), 1)

      (for {
        _ <- dataAccess.createWorkflow(descriptor, Nil, descriptor.namespace.workflow.calls, localBackend)
        _ <- dataAccess.updateWorkflowState(descriptor.id, WorkflowRunning)
        _ <- dataAccess.updateStatus(descriptor.id, Seq(a), ExecutionStatus.Running)
        _ <- dataAccess.updateStatus(descriptor.id, Seq(b), ExecutionStatus.NotStarted)
        wma = (new TestWorkflowManagerSystem).workflowManagerActor
        _ = verifyWorkflowState(wma, descriptor.id, WorkflowSucceeded)
        aStatus <- dataAccess.getExecutionStatus(descriptor.id, a)
        _ = aStatus.map(_.executionStatus) shouldEqual Some(ExecutionStatus.Done)
        bStatus <- dataAccess.getExecutionStatus(descriptor.id, b)
        _ = bStatus.map(_.executionStatus) shouldEqual Some(ExecutionStatus.Done)
      } yield ()).futureValue
    }
  }
}