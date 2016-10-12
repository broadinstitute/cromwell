package cromwell.engine.workflow.lifecycle.execution.ejea

import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.core.WorkflowId
import cromwell.core.callcaching.{CallCachingActivity, ReadAndWriteCache}
import cromwell.core.simpleton.WdlValueSimpleton
import cromwell.engine.workflow.lifecycle.execution.EngineJobExecutionActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.HashError
import cromwell.engine.workflow.lifecycle.execution.callcaching.FetchCachedResultsActor.{CachedOutputLookupFailed, CachedOutputLookupSucceeded}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCachingEntryId
import cromwell.engine.workflow.lifecycle.execution.ejea.EngineJobExecutionActorSpec._
import cromwell.engine.workflow.lifecycle.execution.ejea.HasJobSuccessResponse.SuccessfulCallCacheHashes
import wdl4s.values.WdlString

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
        ejhaResponse foreach { ejea ! _ }

        // Send the response from the "Fetch" actor
        val cachedSimpletons = Seq(WdlValueSimpleton("a", WdlString("hullo")), WdlValueSimpleton("b", WdlString("cheerio")))
        val detritusMap = Map("stdout" -> "//somePath")
        val cachedReturnCode = Some(17)
        val sourceCacheDetails = s"${WorkflowId.randomId}:call-someTask:1"
        ejea ! CachedOutputLookupSucceeded(cachedSimpletons, detritusMap, cachedReturnCode, CallCachingEntryId(75), sourceCacheDetails)
        helper.callCacheHitCopyingProbe.expectMsg(CopyOutputsCommand(cachedSimpletons, detritusMap, cachedReturnCode))

        // Check we end up in the right state:
        ejea.stateName should be(BackendIsCopyingCachedOutputs)
        // Check we end up with the right data:
        val expectedData = initialData.copy(hashes = ejhaResponse map {
          case SuccessfulCallCacheHashes => Success(SuccessfulCallCacheHashes)
          case HashError(t) => Failure(t)
          case _ => fail(s"Bad test wiring. We didn't expect $ejhaResponse here")
        })
        ejea.stateData should be(expectedData)
      }

      RestartOrExecuteCommandTuples foreach { case RestartOrExecuteCommandTuple(operationName, restarting, expectedMessage) =>
        s"Correctly receive $name then $operationName the job when it gets a result-fetch failure" in {
          ejea = ejeaInFetchingCachedOutputsFromDatabaseState(restarting)
          // Send the response from the EJHA (if there was one!):
          ejhaResponse foreach {
            ejea ! _
          }

          // Send the response from the "Fetch" actor
          val failureReason = new Exception("You can't handle the truth!")

          ejea ! CachedOutputLookupFailed(CallCachingEntryId(90210), failureReason)
          helper.bjeaProbe.expectMsg(awaitTimeout, expectedMessage)

          // Check we end up in the right state:
          ejea.stateName should be(RunningJob)
          // Check we end up with the right data:
          val expectedData = initialData.copy(hashes = ejhaResponse map {
            case SuccessfulCallCacheHashes => Success(SuccessfulCallCacheHashes)
            case HashError(t) => Failure(t)
            case _ => fail(s"Bad test wiring. We didn't expect $ejhaResponse here")
          })
          ejea.stateData should be(expectedData)
        }
      }
    }
  }

  def initialData = ResponsePendingData(helper.backendJobDescriptor, helper.bjeaProps, None)
  def ejeaInFetchingCachedOutputsFromDatabaseState(restarting: Boolean = false) = helper.buildEJEA(restarting = restarting, callCachingMode = CallCachingActivity(ReadAndWriteCache)).setStateInline(data = initialData)
}
