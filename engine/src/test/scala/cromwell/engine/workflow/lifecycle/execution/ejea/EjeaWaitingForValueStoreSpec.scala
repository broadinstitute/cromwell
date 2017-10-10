package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.ValueStore
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation

class EjeaWaitingForValueStoreSpec extends EngineJobExecutionActorSpec {

  override implicit val stateUnderTest = CheckingJobStore

  "An EJEA in EjeaWaitingForvalueStore state should" should {
    "prepare the job when receiving the output store" in {
      createWaitingForvalueStoreEjea()
      val valueStore = ValueStore.empty
      ejea ! valueStore
      helper.jobPreparationProbe.expectMsg(awaitTimeout, "expecting CallPreparation Start", CallPreparation.Start(valueStore))
      ejea.stateName should be(PreparingJob)
    }
  }

  private def createWaitingForvalueStoreEjea(): Unit = { ejea = helper.buildEJEA(restarting = true).setStateInline(state = WaitingForvalueStore, data = NoData) }
}
