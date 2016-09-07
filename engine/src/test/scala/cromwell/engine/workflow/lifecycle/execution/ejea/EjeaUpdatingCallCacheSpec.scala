package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import EngineJobExecutionActorSpec._
import cromwell.engine.workflow.lifecycle.execution.callcaching.{CallCacheWriteFailure, CallCacheWriteSuccess}
import cromwell.engine.workflow.lifecycle.execution.ejea.HasJobSuccessResponse.SuccessfulCallCacheHashes
import scala.util.Success

class EjeaUpdatingCallCacheSpec extends EngineJobExecutionActorSpec with HasJobSuccessResponse with CanExpectJobStoreWrites {

  override implicit val stateUnderTest = UpdatingCallCache

  "An EJEA in UpdatingCallCache state" should {

    "Update the Job Store after a successful call cache update" in {
      ejea = ejeaInUpdatingCallCacheState
      ejea ! CallCacheWriteSuccess
      expectJobStoreWrite(initialData)
    }

    "Update the Job Store after a failed call cache update" in {
      ejea = ejeaInUpdatingCallCacheState
      ejea ! CallCacheWriteFailure(new Exception("Alas, poor Yorick! I knew him, Horatio: a fellow of infinite jest, of most excellent fancy"))
      expectJobStoreWrite(initialData)
    }
  }

  lazy val initialData = SucceededResponseData(successResponse, Some(Success(SuccessfulCallCacheHashes)))
  def ejeaInUpdatingCallCacheState = helper.buildEJEA().setStateInline(state = UpdatingCallCache, data = initialData)

}
