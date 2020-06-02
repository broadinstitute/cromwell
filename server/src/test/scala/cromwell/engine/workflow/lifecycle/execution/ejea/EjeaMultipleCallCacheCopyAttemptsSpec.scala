package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.core.WorkflowId
import cromwell.core.callcaching._
import cromwell.core.simpleton.WomValueSimpleton
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.NextHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCachingEntryId
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.CacheHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.CachedOutputLookupSucceeded
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import wom.values.WomString

class EjeaMultipleCallCacheCopyAttemptsSpec
  extends EngineJobExecutionActorSpec
//    with HasJobSuccessResponse
    with HasCopyFailureResponses
//    with HasJobFailureResponses
//    with CanExpectJobStoreWrites
//    with CanExpectCacheWrites
    with CanExpectCacheInvalidation
    with CanExpectFetchCachedResults {

  override implicit val stateUnderTest = BackendIsCopyingCachedOutputs

  "An EJEA attempt to call cache copy" should {

    "Keep trying to copy until it finds a successful hit" in {
      ejea = helper
        .buildEJEA(restarting = false, callCachingMode = CallCachingActivity(ReadCache, CallCachingOptions(invalidateBadCacheResults = false)))
        .setStateInline(state = CheckingCallCache, data = ResponsePendingData(
          jobDescriptor = helper.backendJobDescriptor,
          bjeaProps = helper.bjeaProps,
          ejha = Some(helper.ejhaProbe.ref)))

      0.to(5). foreach { expectedCopyAttemptNumber =>

        // Send back a hit and expect the actor to try to fetch it:
        val callCachingEntryId = CallCachingEntryId(75)
        ejea ! CacheHit(callCachingEntryId)
        expectFetchCachedResultsActor(callCachingEntryId)
        eventually {
          ejea.stateName should be(FetchingCachedOutputsFromDatabase)
          ejea.stateData.asInstanceOf[ResponsePendingData].ejeaCacheHit.get.hitNumber should be(expectedCopyAttemptNumber)
        }

        // The Fetch cached results actor responds with the results which need to be copied:
        val cachedSimpletons = Seq(WomValueSimpleton("a", WomString("hullo")), WomValueSimpleton("b", WomString("cheerio")))
        val detritusMap = Map("stdout" -> "//somePath")
        val cachedReturnCode = Some(17)
        val sourceCacheDetails = s"${WorkflowId.randomId()}:call-someTask:1"
        ejea ! CachedOutputLookupSucceeded(cachedSimpletons, detritusMap, cachedReturnCode, callCachingEntryId, sourceCacheDetails)
        helper.callCacheHitCopyingProbe.expectMsg(CopyOutputsCommand(cachedSimpletons, detritusMap, cachedReturnCode))
        eventually {
          ejea.stateName should be(BackendIsCopyingCachedOutputs)
        }

        // Uh oh - the copy fails. This should send us back to the initial "CheckingCallCache" state of the EJEA
        ejea ! copyAttemptFailedResponse(expectedCopyAttemptNumber)
        helper.ejhaProbe.expectMsg(NextHit)
        eventually {
          ejea.stateName should be(CheckingCallCache)
        }

      }
    }



  }



}
