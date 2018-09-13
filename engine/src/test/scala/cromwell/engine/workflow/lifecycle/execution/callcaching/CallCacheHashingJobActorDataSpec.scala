package cromwell.engine.workflow.lifecycle.execution.callcaching

import cats.data.NonEmptyList
import cromwell.backend._
import cromwell.backend.standard.callcaching.StandardFileHashingActor.SingleFileHashRequest
import cromwell.core.TestKitSuite
import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{CallCacheHashingJobActorData, CompleteFileHashingResult, NoFileHashesResult, PartialFileHashingResult}
import org.scalatest.concurrent.Eventually
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpecLike, Matchers}

class CallCacheHashingJobActorDataSpec extends TestKitSuite with FlatSpecLike with BackendSpec with Matchers with Eventually with TableDrivenPropertyChecks {
  behavior of "CallCacheReadingJobActorData"

  val fileHash1 = HashResult(HashKey("key"), HashValue("value"))
  val fileHash2 = HashResult(HashKey("key2"), HashValue("value2"))
  val fileHash3 = HashResult(HashKey("key3"), HashValue("value3"))
  val fileHashRequest1 = SingleFileHashRequest(null, fileHash1.hashKey, null, null)
  val fileHashRequest2 = SingleFileHashRequest(null, fileHash2.hashKey, null, null)
  val fileHashRequest3 = SingleFileHashRequest(null, fileHash3.hashKey, null, null)
  
  val testCases = Table(
    ("dataBefore", "dataAfter", "result"),
    // No fileHashRequestsRemaining
    (
      CallCacheHashingJobActorData(
        10, List.empty, List.empty, None
      ),
      CallCacheHashingJobActorData(
        10, List.empty, List(fileHash1), None
      ),
      Option(NoFileHashesResult)
    ),
    // Last fileHashRequestsRemaining
    (
      CallCacheHashingJobActorData(
        10, List(List(fileHashRequest1)), List.empty, None
      ),
      CallCacheHashingJobActorData(
        10, List.empty, List(fileHash1), None
      ),
      Option(CompleteFileHashingResult(Set(fileHash1), "6A02F950958AEDA3DBBF83FBB306A030"))
    ),
    // Last batch and not last value
    (
      CallCacheHashingJobActorData(
        10, List(List(fileHashRequest1, fileHashRequest2)), List.empty, None
      ),
      CallCacheHashingJobActorData(
        10, List(List(fileHashRequest2)), List(fileHash1), None
      ),
      None
    ),
    // Not last batch but last value of this batch
    (
      CallCacheHashingJobActorData(
        10, List(List(fileHashRequest1), List(fileHashRequest2)), List.empty, None
      ),
      CallCacheHashingJobActorData(
        10, List(List(fileHashRequest2)), List(fileHash1), None
      ),
      Option(PartialFileHashingResult(NonEmptyList.of(fileHash1)))
    ),
    // Not last batch and not last value of this batch
    (
      CallCacheHashingJobActorData(
        10, List(List(fileHashRequest1, fileHashRequest2), List(fileHashRequest3)), List.empty, None
      ),
      CallCacheHashingJobActorData(
        10, List(List(fileHashRequest2), List(fileHashRequest3)), List(fileHash1), None
      ),
      None
    ),
    // Makes sure new hash is added at the front of the list
    (
      CallCacheHashingJobActorData(
        10, List(List(fileHashRequest1, fileHashRequest2), List(fileHashRequest3)), List(fileHash2), None
      ),
      CallCacheHashingJobActorData(
        10, List(List(fileHashRequest2), List(fileHashRequest3)), List(fileHash1, fileHash2), None
      ),
      None
    ),
    // Use the correct batch size
    (
      CallCacheHashingJobActorData(
        1, List(List.empty, List(fileHashRequest3)), List(fileHash2, fileHash3), None
      ),
      CallCacheHashingJobActorData(
        1, List(List(fileHashRequest3)), List(fileHash1, fileHash2, fileHash3), None
      ),
      Option(PartialFileHashingResult(NonEmptyList.of(fileHash1)))
    )
  )
  
  it should "process new file hashes" in {
    forAll(testCases) { case ((oldData, newData, result)) =>
      oldData.withFileHash(fileHash1) shouldBe (newData -> result)
    }
  }
  
  it should "use the batch size to chop up the file requests" in {
    CallCacheHashingJobActorData(2, List(fileHashRequest1, fileHashRequest2, fileHashRequest3), None) shouldBe
      CallCacheHashingJobActorData(
        2, List(List(fileHashRequest1, fileHashRequest2), List(fileHashRequest3)), List.empty, None
      )
  }
}
