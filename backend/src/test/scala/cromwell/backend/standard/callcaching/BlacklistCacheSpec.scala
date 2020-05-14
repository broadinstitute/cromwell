package cromwell.backend.standard.callcaching

import cromwell.core.CacheConfig
import cromwell.services.CallCaching.CallCachingEntryId
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class BlacklistCacheSpec extends FlatSpec with BeforeAndAfterAll with Matchers {
  val hit = CallCachingEntryId(3)
  "The blacklist cache" should "default, blacklist and expire" in {
    val bucket = "foo"
    val cacheConfig = CacheConfig(concurrency = 1, size = Integer.MAX_VALUE, ttl = 1 second)
    val cache = new RootWorkflowBlacklistCache(cacheConfig)
    cache.isBlacklisted(bucket) shouldBe false
    cache.blacklistBucket(bucket)
    cache.isBlacklisted(bucket) shouldBe true

    // Test ttl
    Thread.sleep(5000L)
    cache.isBlacklisted(bucket) shouldBe false
  }
}
