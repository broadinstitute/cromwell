package cromwell.engine.workflow.lifecycle.execution.ejea

import akka.actor.ActorRef
import cromwell.core.callcaching.{CallCachingActivity, ReadCache}
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.NextHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCacheInvalidatedFailure, CallCacheInvalidatedSuccess}
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._

class EjeaInvalidatingCacheEntrySpec extends EngineJobExecutionActorSpec {

  override implicit val stateUnderTest = InvalidatingCacheEntry

  "An EJEA in InvalidatingCacheEntry state" should {

    val invalidationErrorCause = new Exception("blah")
    val invalidateSuccess = CallCacheInvalidatedSuccess(None)
    val invalidateFailure = CallCacheInvalidatedFailure(invalidationErrorCause)

    List(invalidateSuccess, invalidateFailure) foreach { invalidateActorResponse =>
      s"ask the ejha for the next hit when response is $invalidateActorResponse" in {
        ejea = ejeaInvalidatingCacheEntryState(Option(helper.ejhaProbe.ref))
        // Send the response from the invalidate actor
        ejea ! invalidateActorResponse

        helper.bjeaProbe.expectNoMsg(awaitAlmostNothing)
        helper.ejhaProbe.expectMsg(NextHit)
        eventually { ejea.stateName should be(CheckingCallCache) }
        ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None, Option(helper.ejhaProbe.ref), None))
      }

      RestartOrExecuteCommandTuples foreach { case RestartOrExecuteCommandTuple(operationName, restarting, expectedMessage) =>
        s"$operationName a job if there is no ejha when invalidate response is $invalidateActorResponse" in {
          ejea = ejeaInvalidatingCacheEntryState(None, restarting = restarting)
          // Send the response from the invalidate actor
          ejea ! invalidateActorResponse

          helper.bjeaProbe.expectMsg(awaitTimeout, expectedMessage)
          eventually { ejea.stateName should be(RunningJob) }
          ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper. bjeaProps, None, None))
        }
      }
    }
  }

  def standardResponsePendingData(ejha: Option[ActorRef]) = ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None, ejha)
  def ejeaInvalidatingCacheEntryState(ejha: Option[ActorRef], restarting: Boolean = false) = helper.buildEJEA(restarting = restarting, callCachingMode = CallCachingActivity(ReadCache)).setStateInline(state = InvalidatingCacheEntry, data = standardResponsePendingData(ejha))
}
