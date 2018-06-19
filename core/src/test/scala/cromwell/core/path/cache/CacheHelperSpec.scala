package cromwell.core.path.cache

import cats.effect.IO
import com.google.common.cache.CacheBuilder
import org.scalatest.{FlatSpec, Matchers}

import scala.concurrent.duration._
class CacheHelperSpec extends FlatSpec with Matchers {

  behavior of "CacheHelper"

  it should "cache information properly" in {
    val timeout = 5.seconds
    val fiveSecondsCache = CacheBuilder.newBuilder().expireAfterWrite(timeout.length, timeout.unit).build[String, String]()
    val cacheHelper = new CacheHelper[String, String](fiveSecondsCache) {
      var hasBeenRetrieved: Boolean = false
      override protected def retrieve(key: String) = if (!hasBeenRetrieved) {
        hasBeenRetrieved = true
        IO.pure("hello")
      } else throw new Exception("Not cached !")
    }
    
    cacheHelper.getCachedValue("hi").unsafeRunSync() shouldBe "hello"
    Thread.sleep(1.seconds.toMillis)
    cacheHelper.getCachedValue("hi").unsafeRunSync() shouldBe "hello"
    Thread.sleep(5.seconds.toMillis)
    an[Exception] shouldBe thrownBy { cacheHelper.getCachedValue("hi").unsafeRunSync() }
  }

}
