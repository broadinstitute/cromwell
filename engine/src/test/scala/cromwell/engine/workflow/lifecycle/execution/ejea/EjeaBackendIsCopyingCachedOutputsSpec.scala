package cromwell.engine.workflow.lifecycle.execution.ejea

import cats.data.NonEmptyList
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import EngineJobExecutionActorSpec._
import cromwell.core.callcaching.CallCachingMode
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CallCacheHashes, EJHAResponse, HashError}
import cromwell.engine.workflow.lifecycle.execution.callcaching.MetaInfoId
import scala.util.{Failure, Success, Try}
import cromwell.engine.workflow.lifecycle.execution.ejea.HasJobSuccessResponse.SuccessfulCallCacheHashes

class EjeaBackendIsCopyingCachedOutputsSpec extends EngineJobExecutionActorSpec with HasJobSuccessResponse with HasJobFailureResponses with CanExpectJobStoreWrites with CanExpectCacheWrites with CanExpectCacheInvalidation {

  override implicit val stateUnderTest = BackendIsCopyingCachedOutputs

  "An EJEA in BackendIsCopyingCachedOutputs state" should {

    val hashErrorCause = new Exception("blah")
    val hashResultsDataValue = Some(Success(SuccessfulCallCacheHashes))
    val hashErrorDataValue = Some(Failure(hashErrorCause))

    val hashResultsEjhaResponse = Some(SuccessfulCallCacheHashes)
    val hashErrorEjhaResponse = Some(HashError(hashErrorCause))

    case class InitialHashDataAndEjhaResponseCombination(name: String,
                                                         initialHashData: Option[Try[CallCacheHashes]],
                                                         ejhaResponse: Option[EJHAResponse],
                                                         expectedFinalHashData: Option[Try[CallCacheHashes]],
                                                         validForMode: CallCachingMode => Boolean) {
      def hashingCompletedSuccessfully = ejhaResponse == hashResultsEjhaResponse || initialHashData == hashResultsDataValue
    }

    val initialHashDataAndEjhaResponseCombinations = List(
      // The mode will NEVER have "writeToCache" in these scenarios:
      InitialHashDataAndEjhaResponseCombination("no hashes (ever)", None, None, None, mode => !mode.writeToCache),
      InitialHashDataAndEjhaResponseCombination("unsuccessful hashes in initial data", hashErrorDataValue, None, hashErrorDataValue, mode => !mode.writeToCache),
      // The mode will ALWAYS have "writeToCache" in these scenarios:
      InitialHashDataAndEjhaResponseCombination("unsuccessful hashes from EJHA", None, hashErrorEjhaResponse, hashErrorDataValue, mode => mode.writeToCache),
      InitialHashDataAndEjhaResponseCombination("successful hashes from EJHA", None, hashResultsEjhaResponse, hashResultsDataValue, mode => mode.writeToCache),
      InitialHashDataAndEjhaResponseCombination("successful hashes in initial data", hashResultsDataValue, None, hashResultsDataValue, mode => mode.writeToCache)
    )

    CallCachingModes.filter(_.readFromCache) foreach { mode =>
      initialHashDataAndEjhaResponseCombinations filter { _.validForMode(mode) } foreach { case combo @ InitialHashDataAndEjhaResponseCombination(hashComboName, initialHashData, ejhaResponse, finalHashData, _) =>

        val cacheUpdateRequired = combo.hashingCompletedSuccessfully && mode.writeToCache
        val nextStepName = if (cacheUpdateRequired) "Update call cache" else "Update job store"
        s"$nextStepName when call caching is $mode, the EJEA has $hashComboName and then gets a success result" in {
          ejea = ejeaInBackendIsCopyingCachedOutputsState(initialHashData, mode)
          // Send the response from the EJHA (if there was one!):
          ejhaResponse foreach { ejea ! _ }

          // Nothing should happen here:
          helper.jobStoreProbe.expectNoMsg(awaitAlmostNothing)
          helper.callCacheWriteActorCreations should be(NothingYet) // Rely on the await timeout from the previous step to allow time to pass

          // Send the response from the copying actor
          ejea ! successResponse

          if (cacheUpdateRequired) {
            expectCacheWrite(successResponse, finalHashData.get.get)
          } else {
            expectJobStoreWrite(SucceededResponseData(successResponse, finalHashData))
          }
          // A separate check of the final effective call caching mode:
          if (ejhaResponse == hashErrorEjhaResponse || initialHashData == hashErrorDataValue) {
            ejea.underlyingActor.checkEffectiveCallCachingMode should be(mode.withoutWrite)
          } else {
            ejea.underlyingActor.checkEffectiveCallCachingMode should be(mode)
          }
        }

        s"$nextStepName when it gets a success result and it then gets $hashComboName, if call caching is $mode" in {
          ejea = ejeaInBackendIsCopyingCachedOutputsState(initialHashData, mode)
          // Send the response from the copying actor
          ejea ! successResponse

          ejhaResponse foreach { resp =>
            // Nothing should have happened yet:
            helper.jobStoreProbe.expectNoMsg(awaitAlmostNothing)
            helper.callCacheWriteActorCreations should be(NothingYet) // Rely on the await timeout from the previous step to allow time to pass

            // Ok, now send the response from the EJHA (if there was one!):
            ejea ! resp
          }

          if (cacheUpdateRequired) {
            expectCacheWrite(successResponse, finalHashData.get.get)
          } else {
            expectJobStoreWrite(SucceededResponseData(successResponse, finalHashData))
          }
          // A separate check of the final effective call caching mode:
          if (ejhaResponse == hashErrorEjhaResponse || initialHashData == hashErrorDataValue) {
            ejea.underlyingActor.checkEffectiveCallCachingMode should be(mode.withoutWrite)
          } else {
            ejea.underlyingActor.checkEffectiveCallCachingMode should be(mode)
          }
        }

        RestartOrExecuteCommandTuples foreach { case RestartOrExecuteCommandTuple(operationName, restarting, expectedMessage) =>
          s"invalidate a call for caching if backend coping failed when start mode is $operationName, and it was going to receive $hashComboName, if call caching is $mode" in {
            ejea = ejeaInBackendIsCopyingCachedOutputsState(initialHashData, mode, restarting = restarting)
            // Send the response from the copying actor
            ejea ! failureNonRetryableResponse

            expectInvalidateCallCacheActor(cacheId)
            ejea.stateName should be(InvalidatingCacheEntry)
            ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper. bjeaProps, initialHashData, cacheHit))
          }

          s"invalidate a call for caching if backend coping failed when start mode is $operationName (preserving and received hashes) when call caching is $mode, the EJEA has $hashComboName and then gets a success result" in {
            ejea = ejeaInBackendIsCopyingCachedOutputsState(initialHashData, mode, restarting = restarting)
            // Send the response from the EJHA (if there was one!):
            ejhaResponse foreach { ejea ! _ }

            // Nothing should happen here:
            helper.jobStoreProbe.expectNoMsg(awaitAlmostNothing)
            helper.callCacheWriteActorCreations should be(NothingYet) // Rely on the await timeout from the previous step to allow time to pass

            // Send the response from the copying actor
            ejea ! failureNonRetryableResponse

            expectInvalidateCallCacheActor(cacheId)
            ejea.stateName should be(InvalidatingCacheEntry)
            ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper. bjeaProps, finalHashData, cacheHit))
          }
        }
      }
    }
  }

  private val cacheId: MetaInfoId = MetaInfoId(74)
  private val cacheHit = Option(CacheHit(NonEmptyList.of(cacheId)))
  def standardResponsePendingData(hashes: Option[Try[CallCacheHashes]]) = ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, hashes, cacheHit)
  def ejeaInBackendIsCopyingCachedOutputsState(initialHashes: Option[Try[CallCacheHashes]], callCachingMode: CallCachingMode, restarting: Boolean = false) = helper.buildEJEA(restarting = restarting, callCachingMode = callCachingMode).setStateInline(data = standardResponsePendingData(initialHashes))
}
