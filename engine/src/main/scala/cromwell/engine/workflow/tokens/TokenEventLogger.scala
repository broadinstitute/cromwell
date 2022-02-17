package cromwell.engine.workflow.tokens

import akka.event.LoggingAdapter
import com.google.common.cache.CacheBuilder

import scala.concurrent.duration.FiniteDuration

trait TokenEventLogger {
  def flagTokenHog(hogGroup: String): Unit
  def getLimitedGroups: Set[String]

  def outOfTokens(backend: String): Unit
}

case object NullTokenEventLogger extends TokenEventLogger {
  override def flagTokenHog(hogGroup: String): Unit = ()
  override def getLimitedGroups: Set[String] = Set.empty

  override def outOfTokens(backend: String): Unit = ()
}

class CachingTokenEventLogger(log: LoggingAdapter,
                              cacheEntryTTL: FiniteDuration) extends TokenEventLogger {

  log.info("Using CachingTokenEventLogger with an interval of {} to log token queue events", cacheEntryTTL)

  private val groupCache = CacheBuilder.newBuilder()
    .expireAfterWrite(cacheEntryTTL._1, cacheEntryTTL._2)
    .maximumSize(10000)
    .build[String, Object]()

  // Log when rate limiting is first detected, but do not continually log on subsequent clock ticks
  // Potential issue: the `isEmpty` check prevents a cache entry from refreshing its TTL; instead it expires, is dropped, then re-added
  // The interval before re-adding allows new workflows to unexpectedly sneak through.
  override def flagTokenHog(hogGroup: String): Unit = {
    log.info(s"Token Dispenser: The group $hogGroup has reached its job limit and is being rate-limited.")
    groupCache.put(hogGroup, new Object())
  }

  override def getLimitedGroups: Set[String] = {
    import scala.collection.JavaConverters._
    log.info(s"[TokenEventLogger] Responding with limited groups ${groupCache.asMap()}")
    groupCache.asMap().keySet().asScala.toSet
  }

  private val backendCache = CacheBuilder.newBuilder()
    .expireAfterWrite(cacheEntryTTL._1, cacheEntryTTL._2)
    .maximumSize(10000)
    .build[String, Object]()

  override def outOfTokens(backend: String): Unit = {
    if (Option(backendCache.getIfPresent(backend)).isEmpty) {
      log.info(s"Token Dispenser: The backend $backend is starting too many jobs. New jobs are being limited.")
      backendCache.put(backend, new Object())
    }
  }

}
