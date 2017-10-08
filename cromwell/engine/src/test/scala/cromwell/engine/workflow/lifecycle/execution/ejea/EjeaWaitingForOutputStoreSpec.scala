package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.OutputStore
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation

class EjeaWaitingForOutputStoreSpec extends EngineJobExecutionActorSpec {

  override implicit val stateUnderTest = CheckingJobStore

  "An EJEA in EjeaWaitingForOutputStore state should" should {
    "prepare the job when receiving the output store" in {
      createWaitingForOutputStoreEjea()
      val outputStore = OutputStore.empty
      ejea ! outputStore
      helper.jobPreparationProbe.expectMsg(awaitTimeout, "expecting CallPreparation Start", CallPreparation.Start(outputStore))
      ejea.stateName should be(PreparingJob)
    }
  }

  private def createWaitingForOutputStoreEjea(): Unit = { ejea = helper.buildEJEA(restarting = true).setStateInline(state = WaitingForOutputStore, data = NoData) }
}
