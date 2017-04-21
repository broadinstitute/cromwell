package cromwell.engine.workflow.lifecycle.execution.callcaching

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
  val fileHashRequest1 = mock[SingleFileHashRequest]
  fileHashRequest1.hashKey returns fileHash1.hashKey
  val fileHashRequest2 = mock[SingleFileHashRequest]
  fileHashRequest2.hashKey returns fileHash2.hashKey
  val fileHashRequest3 = mock[SingleFileHashRequest]
  fileHashRequest3.hashKey returns fileHash3.hashKey
  
  val testCases = Table(
    ("dataBefore", "dataAfter", "result"),
    // No fileHashRequestsRemaining
    (
      CallCacheHashingJobActorData(
        List.empty, List.empty, None
      ),
      CallCacheHashingJobActorData(
        List.empty, List(fileHash1), None
      ),
      Option(NoFileHashesResult)
    ),
    // Last fileHashRequestsRemaining
    (
      CallCacheHashingJobActorData(
        List(List(fileHashRequest1)), List.empty, None
      ),
      CallCacheHashingJobActorData(
        List.empty, List(fileHash1), None
      ),
      Option(CompleteFileHashingResult(Set(fileHash1), "2063C1608D6E0BAF80249C42E2BE5804"))
    ),
    // Last batch and not last value
    (
      CallCacheHashingJobActorData(
        List(List(fileHashRequest1, fileHashRequest2)), List.empty, None
      ),
      CallCacheHashingJobActorData(
        List(List(fileHashRequest2)), List(fileHash1), None
      ),
      None
    ),
    // Not last batch but last value of this batch
    (
      CallCacheHashingJobActorData(
        List(List(fileHashRequest1), List(fileHashRequest2)), List.empty, None
      ),
      CallCacheHashingJobActorData(
        List(List(fileHashRequest2)), List(fileHash1), None
      ),
      Option(PartialFileHashingResult(Set(fileHash1)))
    ),
    // Not last batch and not last value of this batch
    (
      CallCacheHashingJobActorData(
        List(List(fileHashRequest1, fileHashRequest2), List(fileHashRequest3)), List.empty, None
      ),
      CallCacheHashingJobActorData(
        List(List(fileHashRequest2), List(fileHashRequest3)), List(fileHash1), None
      ),
      None
    ),
    // Makes sure new hash is added at the front of the list
    (
      CallCacheHashingJobActorData(
        List(List(fileHashRequest1, fileHashRequest2), List(fileHashRequest3)), List(fileHash2), None
      ),
      CallCacheHashingJobActorData(
        List(List(fileHashRequest2), List(fileHashRequest3)), List(fileHash1, fileHash2), None
      ),
      None
    )
  )
  
  it should "process new file hashes" in {
    forAll(testCases) { case ((oldData, newData, result)) =>
      oldData.withFileHash(fileHash1) shouldBe (newData -> result)
    }
  }
}
