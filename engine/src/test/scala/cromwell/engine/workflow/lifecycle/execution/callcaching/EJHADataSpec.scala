package cromwell.engine.workflow.lifecycle.execution.callcaching

import cromwell.core.callcaching._
import cromwell.engine.workflow.lifecycle.execution.callcaching.CallCacheReadActor.CacheResultMatchesForHashes
import org.scalatest.{FlatSpec, Matchers}

class EJHADataSpec extends FlatSpec with Matchers {

  behavior of "EJHA data"

  val hashKey1 = HashKey("Peter Piper picked a peck of pickled peppers")
  val hashKey2 = HashKey("I saw Susie sitting in a shoe shine shop: Where she sits she shines, and where she shines she sits")
  val hashKey3 = HashKey("How many boards could the Mongols hoard if the Mongol hordes got bored?")
  val hashKey4 = HashKey("The sixth sick Sikh's sixth sheep is sick.")
  val hashKey5 = HashKey("The Doge did what a Doge does, when a Doge does his duty to a Duke, that is. When the Doge did his duty and the Duke didn't, that's when the Duchess did the dirt to the Duke with the Doge. There they were in the dark: The Duke with his dagger, the Doge with his dart and the Duchess with her dirk. The Duchess dug at the Duke just when the Duke dove at the Doge. Now the Duke ducked, the Doge dodged, and the Duchess didn't. So the Duke got the Duchess, the Duchess got the Doge, and the Doge got the Duke.")

  val allHashKeys = Set(hashKey1, hashKey2, hashKey3, hashKey4, hashKey5)

  it should "create lists appropriately in the apply method" in {
    val readWriteData = EJHAData(allHashKeys, CallCachingActivity(ReadAndWriteCache))
    readWriteData.remainingCacheChecks should be(allHashKeys)
    readWriteData.remainingHashesNeeded should be(allHashKeys)

    val readOnlyData = EJHAData(allHashKeys, CallCachingActivity(ReadCache))
    readOnlyData.remainingCacheChecks should be(allHashKeys)
    readOnlyData.remainingHashesNeeded should be(Set.empty[HashKey])

    val writeOnlyData = EJHAData(allHashKeys, CallCachingActivity(WriteCache))
    writeOnlyData.remainingCacheChecks should be(Set.empty[HashKey])
    writeOnlyData.remainingHashesNeeded should be(allHashKeys)
  }


  it should "accumulate new hashes" in {
    val data = EJHAData(allHashKeys, CallCachingActivity(WriteCache))
    data.hashesKnown should be(Set.empty)
    data.remainingHashesNeeded should be(allHashKeys)
    data.allHashesKnown should be(false)

    val mostHashKeys = Set(hashKey1, hashKey2, hashKey4, hashKey5)
    val hashResults = mostHashKeys map { x => Set(makeHashResult(x)) }
    val newData = hashResults.foldLeft(data)( (d, h) => d.withNewKnownHashes(h) )
    newData.hashesKnown.map(_.hashKey) should be(mostHashKeys)
    newData.remainingHashesNeeded should be(Set(hashKey3))
    newData.allHashesKnown should be(false)

    val newerData = newData.withNewKnownHashes(Set(makeHashResult(hashKey3)))
    newerData.hashesKnown.map(_.hashKey) should be(allHashKeys)
    newerData.remainingHashesNeeded should be(Set.empty)
    newerData.allHashesKnown should be(true)
  }

  it should "intersect new cache meta info result IDs for cache hits" in {
    val data = EJHAData(allHashKeys, CallCachingActivity(ReadCache))
    data.possibleCacheResults should be(None)
    data.allCacheResultsIntersected should be(false)
    data.isDefinitelyCacheHit should be(false)
    data.isDefinitelyCacheMiss should be(false)

    // To save you time I'll just tell you: the intersection of all these sets is Set(5)
    val cacheLookupResults: List[CacheResultMatchesForHashes] = List(
      CacheResultMatchesForHashes(Set(makeHashResult(hashKey1)), Set(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) map CallCachingEntryId),
      CacheResultMatchesForHashes(Set(makeHashResult(hashKey2)), Set(1, 2, 3, 4, 5, 6) map CallCachingEntryId),
      CacheResultMatchesForHashes(Set(makeHashResult(hashKey3)), Set(1, 2, 3, 5, 7, 8, 9, 10) map CallCachingEntryId),
      CacheResultMatchesForHashes(Set(makeHashResult(hashKey4)), Set(4, 5, 6, 7, 8, 9, 10) map CallCachingEntryId),
      CacheResultMatchesForHashes(Set(makeHashResult(hashKey5)), Set(1, 2, 5, 6, 7, 10) map CallCachingEntryId))
    val newData = cacheLookupResults.foldLeft(data)( (d, c) => d.intersectCacheResults(c) )
    newData.possibleCacheResults match{
      case Some(set) => set should be(Set(CallCachingEntryId(5)))
      case None => fail("There should be a cache result set")
    }
    newData.allCacheResultsIntersected should be(true)
    newData.isDefinitelyCacheHit should be(true)
    newData.isDefinitelyCacheMiss should be(false)
  }

  it should "intersect new cache meta info result IDs for cache misses" in {
    val data = EJHAData(allHashKeys, CallCachingActivity(ReadCache))

    // To save you time I'll just tell you: the intersection of all these sets is empty Set()
    val cacheLookupResults: List[CacheResultMatchesForHashes] = List(
      CacheResultMatchesForHashes(Set(makeHashResult(hashKey1)), Set(1, 2, 3, 4, 5, 6) map CallCachingEntryId),
      CacheResultMatchesForHashes(Set(makeHashResult(hashKey2)), Set(1, 2, 3, 7, 8, 9) map CallCachingEntryId),
      CacheResultMatchesForHashes(Set(makeHashResult(hashKey3)), Set(5, 7, 8, 9, 10) map CallCachingEntryId))
    val newData = cacheLookupResults.foldLeft(data)( (d, c) => d.intersectCacheResults(c) )
    newData.possibleCacheResults match{
      case Some(set) => set should be(Set.empty)
      case None => fail("There should be a cache result set")
    }
    newData.allCacheResultsIntersected should be(false)
    newData.isDefinitelyCacheHit should be(false)
    newData.isDefinitelyCacheMiss should be(true)
  }

  private def makeHashResult(hashKey: HashKey): HashResult = HashResult(hashKey, HashValue("whatever"))
}
