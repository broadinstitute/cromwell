package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.core.WorkflowId
import cromwell.core.callcaching.{CallCachingActivity, ReadAndWriteCache}
import cromwell.core.simpleton.WomValueSimpleton
import cromwell.engine.workflow.lifecycle.execution.job.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, HashError}
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedOutputLookupFailed, CachedOutputLookupSucceeded}
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import cromwell.engine.workflow.lifecycle.execution.ejea.HasJobSuccessResponse.SuccessfulCallCacheHashes
import wom.values.WomString
import common.assertion.CaseClassAssertions._
import cromwell.services.CallCaching.CallCachingEntryId

import scala.util.{Failure, Success}

class EjeaFetchingCachedOutputsFromDatabaseSpec extends EngineJobExecutionActorSpec with HasJobSuccessResponse {

  implicit override def stateUnderTest = FetchingCachedOutputsFromDatabase

  "An EJEA in FetchingCachedOutputsFromDatabase state" should {

    val possibleEjhaResponses = List(
      ("no hashes", None),
      ("hash results", Some(SuccessfulCallCacheHashes)),
      ("hash error", Some(HashError(new Exception("blah")))))

    possibleEjhaResponses foreach { case (name, ejhaResponse) =>
      s"Correctly receive $name then begin the backend-specific output copying actor when it gets a result-fetch success" in {
        ejea = ejeaInFetchingCachedOutputsFromDatabaseState()
        // Send the response from the EJHA (if there was one!):
        ejhaResponse foreach {
          ejea ! _
        }

        // Send the response from the "Fetch" actor
        val cachedSimpletons = Seq(WomValueSimpleton("a", WomString("hullo")), WomValueSimpleton("b", WomString("cheerio")))
        val detritusMap = Map("stdout" -> "//somePath")
        val cachedReturnCode = Some(17)
        val sourceCacheDetails = s"${WorkflowId.randomId()}:call-someTask:1"
        ejea ! CachedOutputLookupSucceeded(cachedSimpletons, detritusMap, cachedReturnCode, callCachingEntryId, sourceCacheDetails)
        helper.callCacheHitCopyingProbe.expectMsg(CopyOutputsCommand(cachedSimpletons, detritusMap, callCachingEntryId, cachedReturnCode))

        // Check we end up in the right state:
        ejea.stateName should be(BackendIsCopyingCachedOutputs)
        // Check we end up with the right data:
        val expectedData = initialData.copy(
          hashes = ejhaResponse map {
            case SuccessfulCallCacheHashes => Success(SuccessfulCallCacheHashes)
            case HashError(t) => Failure(t)
            case _ => fail(s"Bad test wiring. We didn't expect $ejhaResponse here")
          },
          // NB: the "callCachingEntryId" is verified when the EJEA gets the fetch result (and so not changed).
          // Similarly, hit number is not incremented by the FetchingCachedOutputs state (because it was already done):
          // But, we should check that the details are correctly updated:
          ejeaCacheHit = initialData.ejeaCacheHit.map(oldHit => oldHit.copy(details = Some(sourceCacheDetails)))
        )

        ejea.stateData.asInstanceOf[ResponsePendingData] shouldEqualFieldwise expectedData
      }

      RestartOrExecuteCommandTuples foreach { case RestartOrExecuteCommandTuple(operationName, restarting, expectedMessage) =>
        // Send the response from the "Fetch" actor
        val failureReason = new Exception("You can't handle the truth!")
        val lookupFailedMsg = CachedOutputLookupFailed(CallCachingEntryId(90210), failureReason)

        s"Correctly receive $name then $operationName the job when it gets a ${lookupFailedMsg.getClass.getSimpleName} result-fetch failure" in {
          ejea = ejeaInFetchingCachedOutputsFromDatabaseState(restarting)
          // Send the response from the EJHA (if there was one!):
          ejhaResponse foreach {
            ejea ! _
          }

          ejea ! lookupFailedMsg

          helper.bjeaProbe.expectMsg(awaitTimeout, expectedMessage)

          // Check we end up in the right state:
          ejea.stateName should be(RunningJob)
          // Check we end up with the right data:
          val expectedData = initialData.copy(
            hashes = ejhaResponse map {
              case SuccessfulCallCacheHashes => Success(SuccessfulCallCacheHashes)
              case HashError(t) => Failure(t)
              case _ => fail(s"Bad test wiring. We didn't expect $ejhaResponse here")
            },
            backendJobActor = Option(helper.bjeaProbe.ref)
          )
          ejea.stateData should be(expectedData)
        }
      }
    }
  }

  val callCachingEntryId = CallCachingEntryId(75)
  def initialData = ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None, None, Some(EJEACacheHit(CacheHit(callCachingEntryId), 2, None)), None)
  def ejeaInFetchingCachedOutputsFromDatabaseState(restarting: Boolean = false) = helper.buildEJEA(restarting = restarting, callCachingMode = CallCachingActivity(ReadAndWriteCache)).setStateInline(data = initialData)
}
