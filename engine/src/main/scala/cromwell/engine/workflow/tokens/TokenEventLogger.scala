package cromwell.engine.workflow.tokens

import akka.event.LoggingAdapter
import com.google.common.cache.CacheBuilder

import scala.concurrent.duration.FiniteDuration

trait TokenEventLogger {
  def flagTokenHog(hogGroup: String): Unit
  def getLimitedGroups: Set[String]
  def logLimitedGroups(): Unit

  def outOfTokens(backend: String): Unit
  def logOutOfTokens(): Unit
}

case object NullTokenEventLogger extends TokenEventLogger {
  override def flagTokenHog(hogGroup: String): Unit = ()
  override def getLimitedGroups: Set[String] = Set.empty
  override def logLimitedGroups(): Unit = ()

  override def outOfTokens(backend: String): Unit = ()
  override def logOutOfTokens(): Unit = ()
}

class CachingTokenEventLogger(log: LoggingAdapter,
                              cacheEntryTTL: FiniteDuration) extends TokenEventLogger {

  log.info("Using CachingTokenEventLogger with an interval of {} to log token queue events", cacheEntryTTL)

  private val groupCache = CacheBuilder.newBuilder()
    .expireAfterWrite(cacheEntryTTL._1, cacheEntryTTL._2)
    .maximumSize(10000)
    .build[String, Object]()

  override def flagTokenHog(hogGroup: String): Unit = {
    groupCache.put(hogGroup, new Object())
  }

  override def getLimitedGroups: Set[String] = {
    import scala.collection.JavaConverters._
    groupCache.asMap().keySet().asScala.toSet
  }

  override def logLimitedGroups(): Unit = {
    getLimitedGroups foreach { group =>
      log.info(s"Token Dispenser: The group $group has reached its job limit and is being rate-limited.")
    }
  }

  private val backendCache = CacheBuilder.newBuilder()
    .expireAfterWrite(cacheEntryTTL._1, cacheEntryTTL._2)
    .maximumSize(10000)
    .build[String, Object]()

  override def outOfTokens(backend: String): Unit = {
    backendCache.put(backend, new Object())
  }

  override def logOutOfTokens(): Unit = {
    import scala.collection.JavaConverters._
    val backends: Set[String] = backendCache.asMap().keySet().asScala.toSet
    backends foreach { backend =>
      log.info(s"Token Dispenser: The backend $backend is starting too many jobs. New jobs are being limited.")
    }
  }

}
