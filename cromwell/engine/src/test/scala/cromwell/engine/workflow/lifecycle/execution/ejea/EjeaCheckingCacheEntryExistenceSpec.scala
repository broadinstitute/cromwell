package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.core.callcaching.CallCachingOff
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.RequestOutputStore
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec.EnhancedTestEJEA

class EjeaCheckingCacheEntryExistenceSpec extends EngineJobExecutionActorSpec {

  override implicit val stateUnderTest = CheckingJobStore

  "An EJEA in EjeaCheckingCacheEntryExistence state should" should {
    "disable call caching and prepare job if a cache entry already exists for this job" in {
      createCheckingCacheEntryExistenceEjea()
      
      ejea ! HasCallCacheEntry(CallCacheEntryForCall(helper.workflowId, helper.jobDescriptorKey))
      helper.replyToProbe.expectMsg(RequestOutputStore)
      ejea.stateName should be(WaitingForOutputStore)
      
      ejea.underlyingActor.effectiveCallCachingMode shouldBe CallCachingOff
    }

    "prepare job if no cache entry already exists" in {
      createCheckingCacheEntryExistenceEjea()

      ejea ! NoCallCacheEntry(CallCacheEntryForCall(helper.workflowId, helper.jobDescriptorKey))
      helper.replyToProbe.expectMsg(RequestOutputStore)
      ejea.stateName should be(WaitingForOutputStore)
    }

    "prepare job if cache entry existence lookup fails" in {
      createCheckingCacheEntryExistenceEjea()

      ejea ! CacheResultLookupFailure(new Exception("[TEST] Failed to lookup cache entry existence"))
      helper.replyToProbe.expectMsg(RequestOutputStore)
      ejea.stateName should be(WaitingForOutputStore)
    }
  }

  private def createCheckingCacheEntryExistenceEjea(): Unit = { ejea = helper.buildEJEA(restarting = true).setStateInline(state = CheckingCacheEntryExistence, data = NoData) }
}
