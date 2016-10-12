package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.backend.BackendJobExecutionActor.JobFailedNonRetryableResponse
import cromwell.core.callcaching.CallCachingMode
import cromwell.engine.workflow.lifecycle.execution.CallPreparationActor.{BackendJobPreparationSucceeded, CallPreparationFailed}
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import org.scalatest.concurrent.Eventually

class EjeaPreparingJobSpec extends EngineJobExecutionActorSpec with CanExpectHashingInitialization with Eventually {

  override implicit val stateUnderTest = PreparingJob

  "An EJEA in PreparingJob state" should {

    CallCachingModes foreach { mode =>
      if (mode.readFromCache) {
        s"Start checking for a cache hit when job preparation succeeds ($mode)" in {
          ejea = ejeaInPreparingState(mode)
          ejea ! jobPrepSuccessResponse
          expectHashingActorInitialization(mode)
          ejea.stateName should be(CheckingCallCache)
          ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None))
        }
      } else {
        RestartOrExecuteCommandTuples foreach { case RestartOrExecuteCommandTuple(operationName, restarting, expectedMessage) =>
          s"Send BJEA '$operationName' when job preparation succeeds ($mode)" in {
            ejea = ejeaInPreparingState(mode = mode, restarting = restarting)
            ejea ! jobPrepSuccessResponse
            helper.bjeaProbe.expectMsg(awaitTimeout, "job preparation", expectedMessage)
            ejea.stateName should be(RunningJob)
            ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None))
          }
        }
      }

      s"Not proceed if Job Preparation fails ($mode)" in {
        val prepActorResponse = CallPreparationFailed(helper.jobDescriptorKey, new Exception("The goggles! They do nothing!"))
        val prepFailedEjeaResponse = JobFailedNonRetryableResponse(helper.jobDescriptorKey, prepActorResponse.throwable, None)
        ejea = ejeaInPreparingState(mode)
        ejea ! prepActorResponse
        helper.replyToProbe.expectMsg(prepFailedEjeaResponse)
        helper.deathwatch.expectTerminated(ejea, awaitTimeout)
      }
    }
  }

  def jobPrepSuccessResponse = BackendJobPreparationSucceeded(helper.backendJobDescriptor, helper.bjeaProps)

  def ejeaInPreparingState(mode: CallCachingMode, restarting: Boolean = false) = helper.buildEJEA(restarting = restarting, callCachingMode = mode).setStateInline(state = PreparingJob, data = NoData)

}
