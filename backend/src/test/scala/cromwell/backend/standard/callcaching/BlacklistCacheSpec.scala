package cromwell.backend.standard.callcaching

import cromwell.core.WorkflowId
import org.scalatest.{BeforeAndAfterAll, FlatSpec, Matchers}

import scala.concurrent.duration._
import scala.language.postfixOps

class BlacklistCacheSpec extends FlatSpec with BeforeAndAfterAll with Matchers {
  "The blacklist cache" should "default, blacklist and expire" in {
    val id = WorkflowId.randomId()
    val bucket = "foo"
    val cache = new BlacklistCache(10, 2 seconds)
    cache.isBlacklisted(id, bucket) shouldBe false
    cache.blacklist(id, bucket)
    cache.isBlacklisted(id, bucket) shouldBe true

    // Test ttl
    Thread.sleep(5000L)
    cache.isBlacklisted(id, bucket) shouldBe false
  }
}
