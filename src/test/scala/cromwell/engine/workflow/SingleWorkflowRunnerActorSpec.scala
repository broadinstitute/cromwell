package cromwell.engine.workflow

import com.typesafe.config.ConfigFactory
import cromwell.engine.backend.Backend
import cromwell.{CromwellSpec, CromwellTestkitSpec}
import cromwell.engine.db.DataAccess
import cromwell.util.SampleWdl.ThreeStep

import scala.concurrent.Await
import scala.concurrent.duration.Duration
import scala.language.postfixOps

class SingleWorkflowRunnerActorSpec extends CromwellTestkitSpec("SingleWorkflowRunnerActorSpec") {
  val dataAccess = DataAccess()
  val workflowManagerActor = system.actorOf(WorkflowManagerActor.props(dataAccess, CromwellSpec.BackendInstance))
  val props = SingleWorkflowRunnerActor.props(ThreeStep.wdlSource(), ThreeStep.wdlJson, ThreeStep.rawInputs, workflowManagerActor)

  override def afterAll() {
    super.afterAll()
    Await.result(dataAccess.shutdown(), Duration.Inf)
  }

  "A SingleWorkflowRunnerActor" should {
    "successfully run a workflow" in {
      waitForPattern("workflow finished with status 'Succeeded'.") {
        system.actorOf(props)
      }
    }
  }
}
