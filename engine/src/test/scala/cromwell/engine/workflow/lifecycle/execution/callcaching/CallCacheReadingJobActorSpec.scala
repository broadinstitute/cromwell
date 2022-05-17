package cromwell.engine.workflow.lifecycle.execution.callcaching

import akka.testkit.{TestFSMRef, TestProbe}
import cromwell.core.TestKitSuite
import cromwell.core.callcaching.{HashKey, HashResult, HashValue, HashingFailedMessage}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheHashingJobActor.{CompleteFileHashingResult, InitialHashingResult, NextBatchOfFileHashesRequest, NoFileHashesResult}
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadingJobActor.{CCRJAWithData, WaitingForCacheHitOrMiss, _}
import cromwell.engine.workflow.lifecycle.execution.callcaching.EngineJobHashingActor.{CacheHit, CacheMiss, HashError}
import cromwell.services.CallCaching.CallCachingEntryId
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

class CallCacheReadingJobActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Eventually {
  behavior of "CallCacheReadingJobActor"

  it should "try to match initial hashes against DB" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None))
    actorUnderTest.stateName shouldBe WaitingForInitialHash

    // The actual hashes don't matter here, we only care about the aggregated hash
    val aggregatedInitialhash: String = "AggregatedInitialHash"
    callCacheHashingActor.send(actorUnderTest, InitialHashingResult(Set.empty, aggregatedInitialhash))
    callCacheReadProbe.expectMsg(HasMatchingInitialHashLookup(aggregatedInitialhash))
    eventually {
      actorUnderTest.stateName shouldBe WaitingForHashCheck
      actorUnderTest.stateData shouldBe CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialhash, None, Set.empty)
    }
  }

  it should "ask for file hashes if it found matching entries for initial aggregated hash" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None))
    actorUnderTest.setState(WaitingForHashCheck, CCRJAWithData(callCacheHashingActor.ref, "AggregatedInitialHash", None, Set.empty))

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

    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None), parent.ref)
    parent.watch(actorUnderTest)

    actorUnderTest.setState(WaitingForHashCheck, CCRJAWithData(callCacheHashingActor.ref, "AggregatedInitialHash", None, Set.empty))

    callCacheReadProbe.send(actorUnderTest, NoMatchingEntries)
    parent.expectMsg(CacheMiss)
    parent.expectTerminated(actorUnderTest)
  }

  it should "ask for matching cache entries for both aggregated hashes when got both" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()

    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None), TestProbe().ref)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    val aggregatedFileHash: String = "AggregatedFileHash"
    actorUnderTest.setState(WaitingForFileHashes, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, Set.empty))

    val fileHashes = Set(HashResult(HashKey("f1"), HashValue("h1")), HashResult(HashKey("f2"), HashValue("h2")))
    callCacheHashingActor.send(actorUnderTest, CompleteFileHashingResult(fileHashes, aggregatedFileHash))
    callCacheReadProbe.expectMsg(CacheLookupRequest(AggregatedCallHashes(aggregatedInitialHash, aggregatedFileHash), Set.empty, prefixesHint = None))

    eventually {
      actorUnderTest.stateName shouldBe WaitingForCacheHitOrMiss
      actorUnderTest.stateData shouldBe CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, Some(aggregatedFileHash), Set.empty)
    }
  }

  it should "ask for matching cache entries for initial hashes when there is no file input" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()

    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None), TestProbe().ref)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForFileHashes, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, Set.empty))

    callCacheHashingActor.send(actorUnderTest, NoFileHashesResult)
    callCacheReadProbe.expectMsg(CacheLookupRequest(AggregatedCallHashes(aggregatedInitialHash, None), Set.empty, prefixesHint = None))

    eventually {
      actorUnderTest.stateName shouldBe WaitingForCacheHitOrMiss
      actorUnderTest.stateData shouldBe CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, Set.empty)
    }
  }

  it should "reply with next hit when cache hit is successful" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val parent = TestProbe()

    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None), parent.ref)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, Set.empty))

    val id: CallCachingEntryId = CallCachingEntryId(8)
    callCacheReadProbe.send(actorUnderTest, CacheLookupNextHit(id))
    parent.expectMsg(CacheHit(id))

    eventually {
      actorUnderTest.stateName shouldBe WaitingForCacheHitOrMiss
      actorUnderTest.stateData shouldBe CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, Set(id))
    }
  }

  it should "reply with cache miss if there's no hit" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val parent = TestProbe()

    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None), parent.ref)
    parent.watch(actorUnderTest)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, Set.empty))

    callCacheReadProbe.send(actorUnderTest, CacheLookupNoHit)
    parent.expectMsg(CacheMiss)

    parent.expectTerminated(actorUnderTest)
  }

  it should "ask callCacheReadActor for next hit when requested (initial hash only)" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()

    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None), TestProbe().ref)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    val seenCaches: Set[CallCachingEntryId] = Set(CallCachingEntryId(0))
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, seenCaches))

    actorUnderTest ! NextHit
    callCacheReadProbe.expectMsg(CacheLookupRequest(AggregatedCallHashes(aggregatedInitialHash, None), seenCaches, prefixesHint = None))

    actorUnderTest.stateName shouldBe WaitingForCacheHitOrMiss
  }

  it should "ask callCacheReadActor for next hit when requested (with file hash)" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()

    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None), TestProbe().ref)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    val aggregatedFileHash: String = "AggregatedFileHash"
    val seenCaches: Set[CallCachingEntryId] = Set(CallCachingEntryId(0))
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, Option(aggregatedFileHash), seenCaches))

    actorUnderTest ! NextHit
    callCacheReadProbe.expectMsg(CacheLookupRequest(AggregatedCallHashes(aggregatedInitialHash, Option(aggregatedFileHash)), seenCaches, prefixesHint = None))

    actorUnderTest.stateName shouldBe WaitingForCacheHitOrMiss
  }

  it should "reply with cache miss if there's a hash failure" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val parent = TestProbe()

    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None), parent.ref)
    parent.watch(actorUnderTest)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, Set.empty))

    callCacheHashingActor.send(actorUnderTest, HashingFailedMessage("file", new Exception("Hashing failed")))
    parent.expectMsg(CacheMiss)

    parent.expectTerminated(actorUnderTest)
  }

  it should "reply with cache miss if there's a lookup failure" in {
    val callCacheReadProbe = TestProbe()
    val callCacheHashingActor = TestProbe()
    val parent = TestProbe()

    val actorUnderTest = TestFSMRef(new CallCacheReadingJobActor(callCacheReadProbe.ref, prefixesHint = None), parent.ref)
    parent.watch(actorUnderTest)

    val aggregatedInitialHash: String = "AggregatedInitialHash"
    actorUnderTest.setState(WaitingForCacheHitOrMiss, CCRJAWithData(callCacheHashingActor.ref, aggregatedInitialHash, None, Set.empty))

    val reason: Exception = new Exception("Lookup failed")
    callCacheHashingActor.send(actorUnderTest, CacheResultLookupFailure(reason))
    parent.expectMsg(HashError(reason))
    parent.expectMsg(CacheMiss)

    parent.expectTerminated(actorUnderTest)
  }
}
