package cromwell.engine.workflow.lifecycle.execution.callcaching

import cats.data.NonEmptyList
import cromwell.backend._
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.core.TestKitSuite
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{
  CallCacheHashingJobActorData,
  CompleteFileHashingResult,
  NoFileHashesResult,
  PartialFileHashingResult
}
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

class CallCacheHashingJobActorDataSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with BackendSpec
    with Matchers
    with Eventually
    with TableDrivenPropertyChecks {
  behavior of "CallCacheReadingJobActorData"

  private val fileHash1 = HashResult(HashKey("key"), HashValue("value"))
  private val fileHash2 = HashResult(HashKey("key2"), HashValue("value2"))
  private val fileHash3 = HashResult(HashKey("key3"), HashValue("value3"))
  private val fileHashRequest1 = SingleFileHashRequest(null, fileHash1.hashKey, null, null)
  private val fileHashRequest2 = SingleFileHashRequest(null, fileHash2.hashKey, null, null)
  private val fileHashRequest3 = SingleFileHashRequest(null, fileHash3.hashKey, null, null)

  private val testCases = Table(
    ("dataBefore", "dataAfter", "result"),
    // No fileHashRequestsRemaining
    (
      CallCacheHashingJobActorData(
        List.empty,
        List.empty,
        None,
        50
      ),
      CallCacheHashingJobActorData(
        List.empty,
        List(fileHash1),
        None,
        50
      ),
      Option(NoFileHashesResult)
    ),
    // Last fileHashRequestsRemaining
    (
      CallCacheHashingJobActorData(
        List(List(fileHashRequest1)),
        List.empty,
        None,
        50
      ),
      CallCacheHashingJobActorData(
        List.empty,
        List(fileHash1),
        None,
        50
      ),
      Option(CompleteFileHashingResult(Set(fileHash1), "6A02F950958AEDA3DBBF83FBB306A030"))
    ),
    // Last batch and not last value
    (
      CallCacheHashingJobActorData(
        List(List(fileHashRequest1, fileHashRequest2)),
        List.empty,
        None,
        50
      ),
      CallCacheHashingJobActorData(
        List(List(fileHashRequest2)),
        List(fileHash1),
        None,
        50
      ),
      None
    ),
    // Not last batch but last value of this batch
    (
      CallCacheHashingJobActorData(
        List(List(fileHashRequest1), List(fileHashRequest2)),
        List.empty,
        None,
        50
      ),
      CallCacheHashingJobActorData(
        List(List(fileHashRequest2)),
        List(fileHash1),
        None,
        50
      ),
      Option(PartialFileHashingResult(NonEmptyList.of(fileHash1)))
    ),
    // Not last batch and not last value of this batch
    (
      CallCacheHashingJobActorData(
        List(List(fileHashRequest1, fileHashRequest2), List(fileHashRequest3)),
        List.empty,
        None,
        50
      ),
      CallCacheHashingJobActorData(
        List(List(fileHashRequest2), List(fileHashRequest3)),
        List(fileHash1),
        None,
        50
      ),
      None
    ),
    // Makes sure new hash is added at the front of the list
    (
      CallCacheHashingJobActorData(
        List(List(fileHashRequest1, fileHashRequest2), List(fileHashRequest3)),
        List(fileHash2),
        None,
        50
      ),
      CallCacheHashingJobActorData(
        List(List(fileHashRequest2), List(fileHashRequest3)),
        List(fileHash1, fileHash2),
        None,
        50
      ),
      None
    )
  )

  it should "process new file hashes" in {
    forAll(testCases) { case (oldData, newData, result) =>
      oldData.withFileHash(fileHash1) shouldBe (newData -> result)
    }
  }
}
