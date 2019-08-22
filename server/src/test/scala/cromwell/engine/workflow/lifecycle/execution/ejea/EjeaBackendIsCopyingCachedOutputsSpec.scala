package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.NextHit
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCachingEntryId
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CallCacheHashes, EJHAResponse, HashError}
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import cromwell.engine.workflow.lifecycle.execution.ejea.HasJobSuccessResponse.SuccessfulCallCacheHashes

import scala.util.control.NoStackTrace
import scala.util.{Failure, Success, Try}

class EjeaBackendIsCopyingCachedOutputsSpec extends EngineJobExecutionActorSpec with HasJobSuccessResponse with HasCopyFailureResponses with HasJobFailureResponses with CanExpectJobStoreWrites with CanExpectCacheWrites with CanExpectCacheInvalidation {

  override implicit val stateUnderTest = BackendIsCopyingCachedOutputs

  "An EJEA in BackendIsCopyingCachedOutputs state" should {

    val hashErrorCause = new Exception("blah") with NoStackTrace
    val hashResultsDataValue = Some(Success(SuccessfulCallCacheHashes))
    val hashErrorDataValue = Some(Failure(hashErrorCause))

    val hashResultsEjhaResponse: Option[CallCacheHashes] = Some(SuccessfulCallCacheHashes)
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
          helper.jobStoreProbe.expectNoMessage(awaitAlmostNothing)

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
            helper.jobStoreProbe.expectNoMessage(awaitAlmostNothing)

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

        if (mode.readFromCache) {
          s"invalidate a call for caching if backend coping failed when it was going to receive $hashComboName, if call caching is $mode" in {
            ejea = ejeaInBackendIsCopyingCachedOutputsState(initialHashData, mode)
            // Send the response from the copying actor
            ejea ! failedToCopyResponse(cacheHitNumber)

            expectInvalidateCallCacheActor(cacheId)
            eventually {
              ejea.stateName should be(InvalidatingCacheEntry)
            }
            ejea.stateData should be(ResponsePendingData(
              helper.backendJobDescriptor,
              helper.bjeaProps,
              initialHashData,
              Option(helper.ejhaProbe.ref),
              cacheHit,
              None,
              1,
            ))
          }

          s"not invalidate a call for caching if backend coping failed when invalidation is disabled, when it was going to receive $hashComboName, if call caching is $mode" in {
            val invalidationDisabledOptions = CallCachingOptions(invalidateBadCacheResults = false, workflowOptionCallCachePrefixes = None)
            val cacheInvalidationDisabledMode = mode match {
              case CallCachingActivity(rw, _) => CallCachingActivity(rw, invalidationDisabledOptions)
              case _ => fail(s"Mode $mode not appropriate for cache invalidation tests")
            }
            ejea = ejeaInBackendIsCopyingCachedOutputsState(initialHashData, cacheInvalidationDisabledMode)
            // Send the response from the copying actor
            ejea ! failedToCopyResponse(cacheHitNumber)

            helper.ejhaProbe.expectMsg(NextHit)
            
            eventually {
              ejea.stateName should be(CheckingCallCache)
            }
            // Make sure we didn't start invalidating anything:
            helper.invalidateCacheActorCreations.hasExactlyOne should be(false)
            ejea.stateData should be(ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, initialHashData, Option(helper.ejhaProbe.ref), cacheHit))
          }

          s"invalidate a call for caching if backend copying failed (preserving and received hashes) when call caching is $mode, the EJEA has $hashComboName and then gets a success result" in {
            ejea = ejeaInBackendIsCopyingCachedOutputsState(initialHashData, mode)
            // Send the response from the EJHA (if there was one!):
            ejhaResponse foreach {
              ejea ! _
            }

            // Nothing should happen here:
            helper.jobStoreProbe.expectNoMessage(awaitAlmostNothing)

            // Send the response from the copying actor
            ejea ! failedToCopyResponse(cacheHitNumber)

            expectInvalidateCallCacheActor(cacheId)
            eventually {
              ejea.stateName should be(InvalidatingCacheEntry)
            }
            ejea.stateData should be(ResponsePendingData(
              helper.backendJobDescriptor,
              helper.bjeaProps,
              finalHashData,
              Option(helper.ejhaProbe.ref),
              cacheHit,
              None,
              1,
            ))
          }
        }
      }
    }
  }

  private val cacheId: CallCachingEntryId = CallCachingEntryId(74)
  private val cacheHitNumber = 3
  private val cacheHit = Option(EJEACacheHit(CacheHit(cacheId), hitNumber = cacheHitNumber, details = None))
  def standardResponsePendingData(hashes: Option[Try[CallCacheHashes]]) = ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, hashes, Option(helper.ejhaProbe.ref), cacheHit)
  def ejeaInBackendIsCopyingCachedOutputsState(initialHashes: Option[Try[CallCacheHashes]], callCachingMode: CallCachingMode, restarting: Boolean = false) = {
    helper
      .buildEJEA(restarting = restarting, callCachingMode = callCachingMode)
      .setStateInline(data = standardResponsePendingData(initialHashes))
  }
}
