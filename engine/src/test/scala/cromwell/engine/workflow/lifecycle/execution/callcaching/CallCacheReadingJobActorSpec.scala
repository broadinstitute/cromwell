package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.core.TestKitSuite
import cromwell.core.callcaching.{HashKey, HashResult, HashValue, HashingFailedMessage}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{CompleteFileHashingResult, InitialHashingResult, NextBatchOfFileHashesRequest, NoFileHashesResult, PartialFileHashingResult}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.{CCRJAWithData, WaitingForCacheHitOrMiss, _}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CacheMiss, HashError}
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}

class CallCacheReadingJobActorSpec extends TestKitSuite with FlatSpecLike with Matchers with Eventually {
  behavior of "CallCacheReadingJobActor"
  
  it should "try to match initial hashes against DB" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHasingActor = TestProbe()
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref))
    actorUnderTest.stateName shouldBe WaitingForInitialHash

    // The actual hashes don't matter here, we only care about the aggregated hash
    val aggregatedInitialhash: String = "AggregatedInitialHash"
    callCacheHasingActor.send(actorUnderTest, InitialHashingResult(Set.empty, aggregatedInitialhash))
    callCacheReadProbe.expectMsg(HasMatchingInitialHashLookup(aggregatedInitialhash))
    eventually {
      actorUnderTest.stateName shouldBe WaitingForHashCheck
      actorUnderTest.stateData shouldBe CCRJAWithData(callCacheHasingActor.ref, aggregatedInitialhash, None, 1)
    }
  }

  it should "ask for file hashes if it found matching entries for initial aggregated hash" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref))
    actorUnderTest.setState(WaitingForHashCheck, CCRJAWithData(callCacheHashingActor.ref, "AggregatedInitialHash", None, 1))

    callCacheReadProbe.send(actorUnderTest, HasMatchingEntries)
    callCacheHashingActor.expectMsg(NextBatchOfFileHashesRequest)
    eventually {
      actorUnderTest.stateName shouldBe WaitingForFileHashes
    }
  }

  it should "cache miss and die if it didn't find any matching entries for initial aggregated hash" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val parent = TestProbe()
    // self acts as the parent
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref), parent.ref)
    parent.watch(actorUnderTest)
    
    actorUnderTest.setState(WaitingForHashCheck, CCRJAWithData(callCacheHashingActor.ref, "AggregatedInitialHash", None, 1))

    callCacheReadProbe.send(actorUnderTest, NoMatchingEntries)
    parent.expectMsg(CacheMiss)
    parent.expectTerminated(actorUnderTest)
  }

  it should "try to match partial file hashes against DB" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    // self acts as the parent
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref), TestProbe().ref)

    actorUnderTest.setState(WaitingForFileHashes)

    val fileHashes = Set(HashResult(HashKey("f1"), HashValue("h1")), HashResult(HashKey("f2"), HashValue("h2")))
    callCacheHashingActor.send(actorUnderTest, PartialFileHashingResult(fileHashes))
    callCacheReadProbe.expectMsg(HasMatchingInputFilesHashLookup(fileHashes))
    
    eventually {
      actorUnderTest.stateName shouldBe WaitingForHashCheck
    }
  }

  it should "ask for matching cache entries for both aggregated hashes when got both" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    // self acts as the parent
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref), TestProbe().ref)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    val aggregatedFileHash: String = "AggregatedFileHash"
    actorUnderTest.setState(WaitingForFileHashes, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, 1))

    val fileHashes = Set(HashResult(HashKey("f1"), HashValue("h1")), HashResult(HashKey("f2"), HashValue("h2")))
    callCacheHashingActor.send(actorUnderTest, CompleteFileHashingResult(fileHashes, aggregatedFileHash))
    callCacheReadProbe.expectMsg(CacheLookupRequest(AggregatedCallHashes(aggregatedInitialHash, aggregatedFileHash), 1))

    eventually {
      actorUnderTest.stateName shouldBe WaitingForCacheHitOrMiss
      actorUnderTest.stateData shouldBe CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, Some(aggregatedFileHash), 1)
    }
  }

  it should "ask for matching cache entries for initial hashes when there is no file input" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    // self acts as the parent
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref), TestProbe().ref)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForFileHashes, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, 1))

    callCacheHashingActor.send(actorUnderTest, NoFileHashesResult)
    callCacheReadProbe.expectMsg(CacheLookupRequest(AggregatedCallHashes(aggregatedInitialHash, None), 1))

    eventually {
      actorUnderTest.stateName shouldBe WaitingForCacheHitOrMiss
      actorUnderTest.stateData shouldBe CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, 1)
    }
  }

  it should "reply with next hit when cache hit is successful" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val parent = TestProbe()
    // self acts as the parent
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref), parent.ref)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, 1))

    val id: CallCachingEntryId = CallCachingEntryId(8)
    callCacheReadProbe.send(actorUnderTest, CacheLookupNextHit(id))
    parent.expectMsg(CacheHit(id))

    eventually {
      actorUnderTest.stateName shouldBe WaitingForCacheHitOrMiss
      actorUnderTest.stateData shouldBe CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, 2)
    }
  }

  it should "reply with cache miss if there's no hit" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val parent = TestProbe()
    // self acts as the parent
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref), parent.ref)
    parent.watch(actorUnderTest)
    
    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, 1))

    callCacheReadProbe.send(actorUnderTest, CacheLookupNoHit)
    parent.expectMsg(CacheMiss)

    parent.expectTerminated(actorUnderTest)
  }

  it should "ask callCacheReadActor for next hit when requested (initial hash only)" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    // self acts as the parent
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref), TestProbe().ref)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, 2))

    actorUnderTest ! NextHit
    callCacheReadProbe.expectMsg(CacheLookupRequest(AggregatedCallHashes(aggregatedInitialHash, None), 2))

    actorUnderTest.stateName shouldBe WaitingForCacheHitOrMiss
  }

  it should "ask callCacheReadActor for next hit when requested (with file hash)" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    // self acts as the parent
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref), TestProbe().ref)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    val aggregatedFileHash: String = "AggregatedFileHash"
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, Option(aggregatedFileHash), 2))

    actorUnderTest ! NextHit
    callCacheReadProbe.expectMsg(CacheLookupRequest(AggregatedCallHashes(aggregatedInitialHash, Option(aggregatedFileHash)), 2))

    actorUnderTest.stateName shouldBe WaitingForCacheHitOrMiss
  }

  it should "reply with cache miss if there's a hash failure" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val parent = TestProbe()
    // self acts as the parent
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref), parent.ref)
    parent.watch(actorUnderTest)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, 1))

    callCacheHashingActor.send(actorUnderTest, HashingFailedMessage("file", new Exception("Hashing failed")))
    parent.expectMsg(CacheMiss)

    parent.expectTerminated(actorUnderTest)
  }

  it should "reply with cache miss if there's a lookup failure" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val parent = TestProbe()
    // self acts as the parent
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref), parent.ref)
    parent.watch(actorUnderTest)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, 1))

    val reason: Exception = new Exception("Lookup failed")
    callCacheHashingActor.send(actorUnderTest, CacheResultLookupFailure(reason))
    parent.expectMsg(HashError(reason))
    parent.expectMsg(CacheMiss)

    parent.expectTerminated(actorUnderTest)
  }
}
