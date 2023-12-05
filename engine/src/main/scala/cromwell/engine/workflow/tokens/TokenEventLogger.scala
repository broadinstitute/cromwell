package cromwell.engine.workflow.tokens

import com.google.common.cache.CacheBuilder

import scala.concurrent.duration.FiniteDuration

trait TokenEventLogger {
  def flagTokenHog(hogGroup: String): Unit
  def tokenExhaustedGroups: Set[String]

  def outOfTokens(backend: String): Unit
  def tokenExhaustedBackends: Set[String]
}

case object NullTokenEventLogger extends TokenEventLogger {
  override def flagTokenHog(hogGroup: String): Unit = ()
  override def tokenExhaustedGroups: Set[String] = Set.empty

  override def outOfTokens(backend: String): Unit = ()
  override def tokenExhaustedBackends: Set[String] = Set.empty
}

class CachingTokenEventLogger(cacheEntryTTL: FiniteDuration) extends TokenEventLogger {

  private val groupCache = CacheBuilder
    .newBuilder()
    .expireAfterWrite(cacheEntryTTL._1, cacheEntryTTL._2)
    .maximumSize(10000)
    .build[String, Object]()

  override def flagTokenHog(hogGroup: String): Unit =
    groupCache.put(hogGroup, new Object())

  override def tokenExhaustedGroups: Set[String] = {
    import scala.jdk.CollectionConverters._
    groupCache.asMap().keySet().asScala.toSet
  }

  private val backendCache = CacheBuilder
    .newBuilder()
    .expireAfterWrite(cacheEntryTTL._1, cacheEntryTTL._2)
    .maximumSize(10000)
    .build[String, Object]()

  override def outOfTokens(backend: String): Unit =
    backendCache.put(backend, new Object())

  override def tokenExhaustedBackends: Set[String] = {
    import scala.jdk.CollectionConverters._
    backendCache.asMap().keySet().asScala.toSet
  }

}
