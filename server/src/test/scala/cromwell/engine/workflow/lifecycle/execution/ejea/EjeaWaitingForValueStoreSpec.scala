package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.engine.workflow.lifecycle.execution.job.preparation.CallPreparation
import cromwell.engine.workflow.lifecycle.execution.stores.ValueStore

class EjeaWaitingForValueStoreSpec extends EngineJobExecutionActorSpec {

  implicit override val stateUnderTest = CheckingJobStore

  "An EJEA in EjeaWaitingForValueStore state should" should {
    "prepare the job when receiving the output store" in {
      createWaitingForValueStoreEjea()
      val valueStore = ValueStore.empty
      ejea ! valueStore
      helper.jobPreparationProbe.expectMsg(awaitTimeout,
                                           "expecting CallPreparation Start",
                                           CallPreparation.Start(valueStore)
      )
      ejea.stateName should be(PreparingJob)
    }
  }

  private def createWaitingForValueStoreEjea(): Unit = ejea =
    helper.buildEJEA(restarting = true).setStateInline(state = WaitingForValueStore, data = NoData)
}
