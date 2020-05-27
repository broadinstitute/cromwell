package cromwell.backend.standard.callcaching

import cromwell.core.CacheConfig
import cromwell.services.CallCaching.CallCachingEntryId
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class BlacklistCacheSpec extends FlatSpec with Matchers with Eventually {
  "The blacklist cache" should "default, blacklist, whitelist and expire" in {
    val hit = CallCachingEntryId(3)
    val bucket = "foo"

    val bucketCacheConfig = CacheConfig(concurrency = 1, size = Integer.MAX_VALUE, ttl = 1 second)
    val hitCacheConfig = CacheConfig(concurrency = 1, size = Integer.MAX_VALUE, ttl = 1 second)
    val cache = new RootWorkflowBlacklistCache(bucketCacheConfig = bucketCacheConfig, hitCacheConfig = hitCacheConfig)
    cache.getBlacklistStatus(bucket) shouldBe UntestedCacheResult
    cache.getBlacklistStatus(hit) shouldBe UntestedCacheResult
    cache.blacklist(bucket)
    cache.blacklist(hit)
    cache.getBlacklistStatus(bucket) shouldBe BadCacheResult
    cache.getBlacklistStatus(hit) shouldBe BadCacheResult

    implicit val patienceConfig = PatienceConfig(timeout = scaled(5.seconds), interval = scaled(1.second))

    // Test ttl
    eventually {
      cache.getBlacklistStatus(bucket) shouldBe UntestedCacheResult
      cache.getBlacklistStatus(hit) shouldBe UntestedCacheResult
    }

    cache.whitelist(bucket)
    cache.whitelist(hit)
    cache.getBlacklistStatus(bucket) shouldBe GoodCacheResult
    cache.getBlacklistStatus(hit) shouldBe GoodCacheResult

    eventually {
      cache.getBlacklistStatus(bucket) shouldBe UntestedCacheResult
      cache.getBlacklistStatus(hit) shouldBe UntestedCacheResult
    }
  }
}
