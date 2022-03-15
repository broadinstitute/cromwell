package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.database.sql.joins.CallCachingJoin
import cromwell.database.sql.tables._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.engine.workflow.tokens.JobTokenDispenserActor.JobTokenRequest
import cromwell.jobstore.JobStoreActor.RegisterJobCompleted
import cromwell.services.metadata.MetadataService.PutMetadataAction

import scala.util.control.NoStackTrace

class EjeaCheckingCacheEntryExistenceSpec extends EngineJobExecutionActorSpec {

  override implicit val stateUnderTest = CheckingJobStore

  "An EJEA in EjeaCheckingCacheEntryExistence state should" should {
    "re-use the results from the cache hit" in {
      createCheckingCacheEntryExistenceEjea()

      ejea ! CallCachingJoin(CallCachingEntry(helper.workflowId.toString, helper.jobFqn, 0, None, None, allowResultReuse = true),
        List(CallCachingHashEntry("runtime attribute: docker", "HASHVALUE")),
        None,
        List.empty,
        List.empty
      )
      helper.serviceRegistryProbe.expectMsgPF(awaitTimeout) {
        case put: PutMetadataAction => put.events.find(_.key.key.endsWith("runtime attribute:docker"))flatMap(_.value.map(_.value)) shouldBe Option("HASHVALUE")
      }
      helper.jobStoreProbe.expectMsgClass(classOf[RegisterJobCompleted])
      ejea.stateName should be(UpdatingJobStore)
    }

    "prepare job if no cache entry already exists" in {
      createCheckingCacheEntryExistenceEjea()

      ejea ! NoCallCacheEntry(CallCacheEntryForCall(helper.workflowId, helper.jobDescriptorKey))
      helper.jobExecutionTokenDispenserProbe.expectMsgClass(max = awaitTimeout, classOf[JobTokenRequest])
      ejea.stateName should be(RequestingExecutionToken)
    }

    "prepare job if cache entry existence lookup fails" in {
      createCheckingCacheEntryExistenceEjea()

      ejea ! CacheResultLookupFailure(new Exception("[TEST] Failed to lookup cache entry existence") with NoStackTrace)
      helper.jobExecutionTokenDispenserProbe.expectMsgClass(max = awaitTimeout, classOf[JobTokenRequest])
      ejea.stateName should be(RequestingExecutionToken)
    }
  }

  private def createCheckingCacheEntryExistenceEjea(): Unit = {
    ejea = helper.buildEJEA().setStateInline(state = CheckingCacheEntryExistence, data = NoData)
  }
}
