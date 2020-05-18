package cromwell.backend.standard.callcaching

import cromwell.core.CacheConfig
import cromwell.services.CallCaching.CallCachingEntryId
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class BlacklistCacheSpec extends FlatSpec with Matchers {
  "The blacklist cache" should "default, blacklist, whitelist and expire" in {
    val hit = CallCachingEntryId(3)
    val bucket = "foo"

    val bucketCacheConfig = CacheConfig(concurrency = 1, size = Integer.MAX_VALUE, ttl = 1 second)
    val hitCacheConfig = CacheConfig(concurrency = 1, size = Integer.MAX_VALUE, ttl = 1 second)
    val cache = new RootWorkflowBlacklistCache(bucketCacheConfig = bucketCacheConfig, hitCacheConfig = hitCacheConfig)
    cache.getBlacklistStatus(bucket) shouldBe Unknown
    cache.getBlacklistStatus(hit) shouldBe Unknown
    cache.blacklist(bucket)
    cache.blacklist(hit)
    cache.getBlacklistStatus(bucket) shouldBe KnownBad
    cache.getBlacklistStatus(hit) shouldBe KnownBad

    // Test ttl
    Thread.sleep(5000L)
    cache.getBlacklistStatus(bucket) shouldBe Unknown
    cache.getBlacklistStatus(hit) shouldBe Unknown

    cache.whitelist(bucket)
    cache.whitelist(hit)
    cache.getBlacklistStatus(bucket) shouldBe KnownGood
    cache.getBlacklistStatus(hit) shouldBe KnownGood

    Thread.sleep(5000L)
    cache.getBlacklistStatus(bucket) shouldBe Unknown
    cache.getBlacklistStatus(hit) shouldBe Unknown
  }
}
