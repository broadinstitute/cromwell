package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor.{CheckingCallCache, FetchingCachedOutputsFromDatabase, ResponsePendingData, RunningJob}
import EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.core.callcaching.{CallCachingActivity, CallCachingOff, ReadCache}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CacheMiss, HashError}
import cromwell.engine.workflow.lifecycle.execution.callcaching.MetaInfoId
import org.scalatest.concurrent.Eventually

class EjeaCheckingCallCacheSpec extends EngineJobExecutionActorSpec with Eventually {

  override implicit val stateUnderTest = CheckingCallCache

  "An EJEA in CheckingCallCache mode" should {
    "Try to fetch the call cache outputs if it gets a CacheHit" in {
      createCheckingCallCacheEjea()
      ejea ! CacheHit(MetaInfoId(75))
      eventually { helper.fetchCachedResultsActorCreations.hasExactlyOne should be(true) }
      helper.fetchCachedResultsActorCreations checkIt {
        case (CacheHit(metainfoId), _) => metainfoId should be(MetaInfoId(75))
        case _ => fail("Incorrect creation of the fetchCachedResultsActor")
      }
      ejea.stateName should be(FetchingCachedOutputsFromDatabase)
    }

    RestartOrExecuteCommandTuples foreach { case RestartOrExecuteCommandTuple(operationName, restarting, expectedMessage) =>
      s"$operationName the job if it receives a cache miss and restarting is $restarting" in {
        createCheckingCallCacheEjea(restarting)
        ejea ! CacheMiss
        helper.bjeaProbe.expectMsg(awaitTimeout, expectedMessage)
        ejea.stateName should be(RunningJob)
      }

      s"Disabling call caching and $operationName the job if it receives a HashError and restarting is $restarting" in {
        createCheckingCallCacheEjea(restarting)
        ejea ! HashError(new Exception("Anticipated exception. Part of test-flow"))
        eventually {
          ejea.underlyingActor.checkEffectiveCallCachingMode should be(CallCachingOff)
        }
        helper.bjeaProbe.expectMsg(awaitTimeout, expectedMessage)
        ejea.stateName should be(RunningJob)
      }
    }
  }

  private def createCheckingCallCacheEjea(restarting: Boolean = false): Unit = {
    ejea = helper.buildEJEA(restarting = restarting, callCachingMode = CallCachingActivity(ReadCache))
    ejea.setStateInline(state = CheckingCallCache, data = ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None))
    ()
  }
}
