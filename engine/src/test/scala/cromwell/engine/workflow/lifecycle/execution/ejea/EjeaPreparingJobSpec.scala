package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.backend.BackendJobDescriptor
import cromwell.backend.BackendJobExecutionActor.JobFailedNonRetryableResponse
import cromwell.core.callcaching.{CallCachingMode, DockerWithHash}
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import cromwell.engine.workflow.lifecycle.execution.preparation.CallPreparation.{BackendJobPreparationSucceeded, CallPreparationFailed}
import org.scalatest.concurrent.Eventually

class EjeaPreparingJobSpec extends EngineJobExecutionActorSpec with CanExpectHashingInitialization with Eventually {

  override implicit val stateUnderTest = PreparingJob

  "An EJEA in PreparingJob state" should {

    CallCachingModes foreach { mode =>
      if (mode.readFromCache) {
        s"Start checking for a cache hit when job preparation succeeds and a docker hash is available ($mode)" in {
          val jobDescriptor = helper.backendJobDescriptor.copy(maybeCallCachingEligible = DockerWithHash("hello"))
          ejea = ejeaInPreparingState(mode)
          ejea ! jobPrepSuccessResponse(jobDescriptor)
          expectHashingActorInitialization(mode, jobDescriptor)
          ejea.stateName should be(CheckingCallCache)
          ejea.stateData should be(ResponsePendingData(jobDescriptor, helper.bjeaProps, None, Option(helper.ejhaProbe.ref)))
        }

        s"Not check for a cache hit when job preparation succeeds and no docker hash is available ($mode)" in {
          ejea = ejeaInPreparingState(mode)
          ejea ! jobPrepSuccessResponse(helper.backendJobDescriptor)
          ejea.stateName should be(RunningJob)
          ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None))
        }
      } else {
        RestartOrExecuteCommandTuples foreach { case RestartOrExecuteCommandTuple(operationName, restarting, expectedMessage) =>
          s"Send BJEA '$operationName' when job preparation succeeds ($mode)" in {
            ejea = ejeaInPreparingState(mode = mode, restarting = restarting)
            ejea ! jobPrepSuccessResponse(helper.backendJobDescriptor)
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

  def jobPrepSuccessResponse(jobDescriptor: BackendJobDescriptor) = BackendJobPreparationSucceeded(jobDescriptor, helper.bjeaProps)

  def ejeaInPreparingState(mode: CallCachingMode, restarting: Boolean = false) = helper.buildEJEA(restarting = restarting, callCachingMode = mode).setStateInline(state = PreparingJob, data = NoData)

}
