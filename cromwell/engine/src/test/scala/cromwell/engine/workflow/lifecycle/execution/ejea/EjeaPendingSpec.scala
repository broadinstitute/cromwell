package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.tokens.JobExecutionTokenDispenserActor.JobExecutionTokenRequest
import org.scalatest.concurrent.Eventually

class EjeaPendingSpec extends EngineJobExecutionActorSpec with CanValidateJobStoreKey with Eventually {

  override implicit val stateUnderTest: EngineJobExecutionActorState = Pending

  "An EJEA in the Pending state" should {

    List(false, true) foreach { restarting =>
      s"wait for the Execute signal then request an execution token (with restarting=$restarting)" in {
        ejea = helper.buildEJEA(restarting = restarting)
        ejea ! Execute

        helper.jobTokenDispenserProbe.expectMsgClass(max = awaitTimeout, classOf[JobExecutionTokenRequest])

        helper.jobPreparationProbe.msgAvailable should be(false)
        helper.jobStoreProbe.msgAvailable should be(false)
        ejea.stateName should be(RequestingExecutionToken)
      }
    }
  }
}
