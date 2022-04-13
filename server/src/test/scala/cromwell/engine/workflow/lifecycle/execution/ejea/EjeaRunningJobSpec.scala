package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.HashError
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec.EnhancedTestEJEA
import cromwell.engine.workflow.lifecycle.execution.ejea.HasJobSuccessResponse.SuccessfulCallCacheHashes
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor.{FailedResponseData, ResponsePendingData, RunningJob, SucceededResponseData}
import org.scalatest.concurrent.Eventually

import scala.util.control.NoStackTrace
import scala.util.{Failure, Success}

class EjeaRunningJobSpec extends EngineJobExecutionActorSpec with Eventually with CanValidateJobStoreKey with CanExpectJobStoreWrites with CanExpectCacheWrites with HasJobSuccessResponse with HasJobFailureResponses {

  override implicit val stateUnderTest = RunningJob

  val hashError = HashError(new Exception("ARGH!!!") with NoStackTrace)

  val failureCases = List(
    ("FailedRetryableResponse", () => failureRetryableResponse, true),
    ("FailedNonRetryableResponse", () => failureNonRetryableResponse, false)
  )

  "A 'RunningJob' EJEA" should {
    CallCachingModes foreach { mode =>

      /* *************************** */
      /* JobAbortedResponse Handling */
      /* *************************** */

      s"Handle receiving a JobAbortedResponse correctly in $mode mode without hashes" in {
        ejea = ejeaInRunningState(mode)
        ejea ! abortedResponse

        helper.replyToProbe.expectMsg(max = awaitTimeout, hint = "parent wants the response", abortedResponse)

        helper.deathwatch.expectTerminated(ejea, awaitTimeout)
        // Make sure nothing was sent to the JobStore or CacheResultSaver in the meanwhile:
        helper.jobStoreProbe.expectNoMessage(awaitAlmostNothing)
        helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
      }

      s"Handle receiving a JobAbortedResponse correctly in $mode mode with successful hashes" in {
        ejea = ejeaInRunningState(mode)
        ejea ! SuccessfulCallCacheHashes
        eventually { ejea.stateData should be(initialData.copy(hashes = Some(Success(SuccessfulCallCacheHashes)))) }
        ejea.stateName should be(RunningJob)

        ejea ! abortedResponse

        helper.replyToProbe.expectMsg(max = awaitTimeout, hint = "parent wants the response", abortedResponse)

        helper.deathwatch.expectTerminated(ejea, awaitTimeout)
        // Make sure nothing was sent to the JobStore or CacheResultSaver in the meanwhile:
        helper.jobStoreProbe.expectNoMessage(awaitAlmostNothing)
        helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
      }

      s"Handle receiving a JobAbortedResponse correctly in $mode mode with failed hashes" in {
        ejea = ejeaInRunningState(mode)
        ejea ! hashError
        eventually { ejea.stateData should be(initialData.copy(hashes = Some(Failure(hashError.reason)))) }
        ejea.stateName should be(RunningJob)

        ejea ! abortedResponse

        helper.replyToProbe.expectMsg(max = awaitTimeout, hint = "parent wants the response", abortedResponse)

        helper.deathwatch.expectTerminated(ejea, awaitTimeout)
        // Make sure nothing was sent to the JobStore or CacheResultSaver in the meanwhile:
        helper.jobStoreProbe.expectNoMessage(awaitAlmostNothing)
        helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
      }

      /* ********************************************* */
      /* JobAbortedResponse/JobFailedResponse Handling */
      /* ********************************************* */

      // Depending on the value of writeToCache the behavior differes a little:
      // If true, we need to wait for a hash response (successful or not) before completing the job, and then write it to the call cache
      // If false, we don't want need to wait for the hashes and even if we have them we don't want to write them to the call cache
      if (mode.writeToCache) {
        /* JobSuccessResponse */

        s"Handle receiving SuccessResponse then CallCacheHashes correctly in $mode mode" in {
          ejea = ejeaInRunningState(mode)
          ejea ! successResponse
          eventually { ejea.stateData should be(SucceededResponseData(successResponse, None)) }
          ejea.stateName should be(RunningJob)
          ejea ! SuccessfulCallCacheHashes
          expectCacheWrite(successResponse, SuccessfulCallCacheHashes)
        }

        s"Handle receiving CallCacheHashes then SuccessResponse correctly in $mode mode" in {
          ejea = ejeaInRunningState(mode)
          ejea ! SuccessfulCallCacheHashes
          eventually { ejea.stateData should be(initialData.copy(hashes = Some(Success(SuccessfulCallCacheHashes)))) }
          ejea.stateName should be(RunningJob)
          ejea ! successResponse
          expectCacheWrite(successResponse, SuccessfulCallCacheHashes)
        }

        s"Handle receiving SuccessResponse then HashError correctly in $mode mode" in {
          ejea = ejeaInRunningState(mode)
          ejea ! successResponse
          eventually { ejea.stateData should be(SucceededResponseData(successResponse, None)) }
          ejea.stateName should be(RunningJob)
          ejea ! hashError
          expectJobStoreWrite(SucceededResponseData(successResponse, Some(Failure(hashError.reason))))
        }

        s"Handle receiving HashError then SuccessResponse correctly in $mode mode" in {
          ejea = ejeaInRunningState(mode)
          ejea ! hashError
          eventually { ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, Some(Failure(hashError.reason)))) }
          ejea.stateName should be(RunningJob)
          ejea ! successResponse
          expectJobStoreWrite(SucceededResponseData(successResponse, Some(Failure(hashError.reason))))
        }

        /* JobFailedResponse */

        failureCases foreach { case (name, responseMaker, retryable) =>
          s"Handle receiving $name then CallCacheHashes correctly in $mode mode" in {
            val failedResponse = responseMaker()
            ejea = ejeaInRunningState(mode)
            ejea ! failedResponse
            eventually { ejea.stateData should be(FailedResponseData(failedResponse, None)) }
            ejea.stateName should be(RunningJob)
            ejea ! SuccessfulCallCacheHashes
            // Don't expect cache write here (the job failed !) but do expect job store write
            expectJobStoreWriteFailed(FailedResponseData(failedResponse, Option(Success(SuccessfulCallCacheHashes))), retryable)
            helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
          }

          s"Handle receiving CallCacheHashes then $name correctly in $mode mode" in {
            val failedResponse = responseMaker()
            ejea = ejeaInRunningState(mode)
            ejea.stateName should be(RunningJob)
            ejea ! SuccessfulCallCacheHashes
            eventually { ejea.stateData should be(initialData.copy(hashes = Some(Success(SuccessfulCallCacheHashes)))) }
            ejea ! failedResponse
            // Don't expect cache write here (the job failed !) but do expect job store write
            expectJobStoreWriteFailed(FailedResponseData(failedResponse, Option(Success(SuccessfulCallCacheHashes))), retryable)
            helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
          }

          s"Handle receiving $name then HashError correctly in $mode mode" in {
            val failedResponse = responseMaker()
            ejea = ejeaInRunningState(mode)
            ejea ! failedResponse
            eventually { ejea.stateData should be(FailedResponseData(failedResponse, None)) }
            ejea.stateName should be(RunningJob)
            ejea ! hashError
            expectJobStoreWriteFailed(FailedResponseData(failedResponse, Some(Failure(hashError.reason))), retryable)
            helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
          }

          s"Handle receiving HashError then $name correctly in $mode mode" in {
            val failedResponse = responseMaker()
            ejea = ejeaInRunningState(mode)
            ejea ! hashError
            eventually { ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, Some(Failure(hashError.reason)))) }
            ejea.stateName should be(RunningJob)
            ejea ! failedResponse
            expectJobStoreWriteFailed(FailedResponseData(failedResponse, Some(Failure(hashError.reason))), retryable)
            helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
          }
        }
        // writeToCache = false
      } else {
        s"Handle receiving a SucceededResponse correctly in $mode mode without hashes" in {
          ejea = ejeaInRunningState(mode)
          ejea ! successResponse
          expectJobStoreWrite(SucceededResponseData(successResponse, None))
          helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
        }

        s"Handle receiving a SucceededResponse correctly in $mode mode with successful hashes" in {
          ejea = ejeaInRunningState(mode)
          ejea ! SuccessfulCallCacheHashes
          eventually { ejea.stateData should be(initialData.copy(hashes = Some(Success(SuccessfulCallCacheHashes)))) }
          ejea.stateName should be(RunningJob)
          ejea ! successResponse
          // Even if we received hashes, writeToCache is false so we go straight to job store and don't write them to the cache
          expectJobStoreWrite(SucceededResponseData(successResponse, Some(Success(SuccessfulCallCacheHashes))))
          helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
        }

        s"Handle receiving a SucceededResponse correctly in $mode mode with failed hashes" in {
          ejea = ejeaInRunningState(mode)
          ejea ! hashError
          eventually { ejea.stateData should be(initialData.copy(hashes = Some(Failure(hashError.reason)))) }
          ejea.stateName should be(RunningJob)
          ejea ! successResponse
          // Even if we received hashes, writeToCache is false so we go straight to job store and don't write them to the cache
          expectJobStoreWrite(SucceededResponseData(successResponse, Some(Failure(hashError.reason))))
          helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
        }

        failureCases foreach { case (name, responseMaker, retryable) =>
          s"Handle receiving a $name correctly in $mode mode without hashes" in {
            val failedResponse = responseMaker()
            ejea = ejeaInRunningState(mode)
            ejea ! failedResponse
            // writeToCache is off, we don't need to wait for hashes
            expectJobStoreWriteFailed(FailedResponseData(failedResponse, None), retryable)
            helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
          }

          s"Handle receiving a $name correctly in $mode mode with successful hashes" in {
            val failedResponse = responseMaker()
            ejea = ejeaInRunningState(mode)
            ejea ! SuccessfulCallCacheHashes
            eventually { ejea.stateData should be(initialData.copy(hashes = Some(Success(SuccessfulCallCacheHashes)))) }
            ejea.stateName should be(RunningJob)
            ejea ! failedResponse
            expectJobStoreWriteFailed(FailedResponseData(failedResponse, Some(Success(SuccessfulCallCacheHashes))), retryable)
            helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
          }

          s"Handle receiving a $name correctly in $mode mode with failed hashes" in {
            val failedResponse = responseMaker()
            ejea = ejeaInRunningState(mode)
            ejea ! hashError
            eventually { ejea.stateData should be(initialData.copy(hashes = Some(Failure(hashError.reason)))) }
            ejea.stateName should be(RunningJob)
            ejea ! failedResponse
            expectJobStoreWriteFailed(FailedResponseData(failedResponse, Some(Failure(hashError.reason))), retryable)
            helper.callCacheWriteActorProbe.expectNoMessage(awaitAlmostNothing)
          }
        }
      }
    }
  }

  def initialData = ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None)
  def ejeaInRunningState(mode: CallCachingMode = CallCachingActivity(ReadAndWriteCache)) = helper.buildEJEA(callCachingMode = mode).setStateInline(state = RunningJob, data = initialData)
}
