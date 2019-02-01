package cromwell.engine.workflow.tokens

import akka.event.LoggingAdapter
import com.google.common.cache.CacheBuilder

import scala.concurrent.duration.FiniteDuration

trait TokenEventLogger {
  def flagTokenHog(hogGroup: String): Unit
  def outOfTokens(backend: String): Unit
}

case object NullTokenEventLogger extends TokenEventLogger {
  override def flagTokenHog(hogGroup: String): Unit = ()
  override def outOfTokens(backend: String): Unit = ()
}

class CachingTokenEventLogger(log: LoggingAdapter,
                              cacheEntryTTL: FiniteDuration) extends TokenEventLogger {

  log.info("Using CachingTokenEventLogger with an interval of {} to log token queue events", cacheEntryTTL)

  private val cache = CacheBuilder.newBuilder()
    .expireAfterWrite(cacheEntryTTL._1, cacheEntryTTL._2)
    .maximumSize(10000)
    .build[String, Object]()


  override def flagTokenHog(hogGroup: String): Unit = {
    if (Option(cache.getIfPresent("HOG_" + hogGroup)).isEmpty) {
      log.info(s"Token Dispenser: The hog group $hogGroup is starting too many jobs. It is being hog limited.")
      cache.put("HOG_" + hogGroup, new Object())
    }
  }


  override def outOfTokens(backend: String): Unit = {
    if (Option(cache.getIfPresent("OOT_" + backend)).isEmpty) {
      log.info(s"Token Dispenser: The backend $backend is starting too many jobs. New jobs are being limited.")
      cache.put("OOT_" + backend, new Object())
    }
  }

}
