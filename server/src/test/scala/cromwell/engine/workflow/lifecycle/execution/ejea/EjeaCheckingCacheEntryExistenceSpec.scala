package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables._
import cromwell.engine.workflow.lifecycle.execution.WorkflowExecutionActor.RequestValueStore
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.jobstore.JobStoreActor.RegisterJobCompleted

class EjeaCheckingCacheEntryExistenceSpec extends EngineJobExecutionActorSpec {

  override implicit val stateUnderTest = CheckingJobStore

  "An EJEA in EjeaCheckingCacheEntryExistence state should" should {
    "re-use the results from the cache hit" in {
      createCheckingCacheEntryExistenceEjea()
      
      ejea ! CallCachingJoin(CallCachingEntry(helper.workflowId.toString, helper.jobFqn, 0, None, None, allowResultReuse = true),
        List.empty,
        None,
        List.empty,
        List.empty
      )
      helper.jobStoreProbe.expectMsgClass(classOf[RegisterJobCompleted])
      ejea.stateName should be(UpdatingJobStore)
    }

    "prepare job if no cache entry already exists" in {
      createCheckingCacheEntryExistenceEjea()

      ejea ! NoCallCacheEntry(CallCacheEntryForCall(helper.workflowId, helper.jobDescriptorKey))
      helper.replyToProbe.expectMsg(RequestValueStore)
      ejea.stateName should be(WaitingForValueStore)
    }

    "prepare job if cache entry existence lookup fails" in {
      createCheckingCacheEntryExistenceEjea()

      ejea ! CacheResultLookupFailure(new Exception("[TEST] Failed to lookup cache entry existence"))
      helper.replyToProbe.expectMsg(RequestValueStore)
      ejea.stateName should be(WaitingForValueStore)
    }
  }

  private def createCheckingCacheEntryExistenceEjea(): Unit = { ejea = helper.buildEJEA(restarting = true).setStateInline(state = CheckingCacheEntryExistence, data = NoData) }
}
