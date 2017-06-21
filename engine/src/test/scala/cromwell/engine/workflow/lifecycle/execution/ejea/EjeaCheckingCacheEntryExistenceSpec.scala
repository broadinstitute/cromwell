package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.core.callcaching.CallCachingOff
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.{CheckingCacheEntryExistence, CheckingJobStore, NoData, PreparingJob}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation

class EjeaCheckingCacheEntryExistenceSpec extends EngineJobExecutionActorSpec {

  override implicit val stateUnderTest = CheckingJobStore

  "An EJEA in EjeaCheckingCacheEntryExistence state should" should {
    "disable call caching and prepare job if a cache entry already exists for this job" in {
      createCheckingCacheEntryExistenceEjea()
      
      ejea ! HasCallCacheEntry(CallCacheEntryForCall(helper.workflowId, helper.jobDescriptorKey))
      helper.jobPreparationProbe.expectMsg(awaitTimeout, "expecting CallPreparation Start", CallPreparation.Start)
      ejea.stateName should be(PreparingJob)
      
      ejea.underlyingActor.effectiveCallCachingMode shouldBe CallCachingOff
    }

    "prepare job if a cache entry already exists for this job" in {
      createCheckingCacheEntryExistenceEjea()

      ejea ! NoCallCacheEntry(CallCacheEntryForCall(helper.workflowId, helper.jobDescriptorKey))

      helper.jobPreparationProbe.expectMsg(awaitTimeout, "expecting CallPreparation Start", CallPreparation.Start)
      ejea.stateName should be(PreparingJob)
    }

    "prepare job if cache entry existence lookup fails" in {
      createCheckingCacheEntryExistenceEjea()

      ejea ! CacheResultLookupFailure(new Exception("[TEST] Failed to lookup cache entry existence"))

      helper.jobPreparationProbe.expectMsg(awaitTimeout, "expecting CallPreparation Start", CallPreparation.Start)
      ejea.stateName should be(PreparingJob)
    }
  }

  private def createCheckingCacheEntryExistenceEjea(): Unit = { ejea = helper.buildEJEA(restarting = true).setStateInline(state = CheckingCacheEntryExistence, data = NoData) }
}
