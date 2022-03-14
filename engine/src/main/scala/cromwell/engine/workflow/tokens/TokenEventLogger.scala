package cromwell.engine.workflow.tokens

import com.google.common.cache.CacheBuilder

import scala.concurrent.duration.FiniteDuration

trait TokenEventLogger {
  def flagTokenHog(hogGroup: String): Unit
  def getLimitedGroups: Set[String]

  def outOfTokens(backend: String): Unit
  def getLimitedBackends: Set[String]
}

case object NullTokenEventLogger extends TokenEventLogger {
  override def flagTokenHog(hogGroup: String): Unit = ()
  override def getLimitedGroups: Set[String] = Set.empty

  override def outOfTokens(backend: String): Unit = ()
  override def getLimitedBackends: Set[String] = Set.empty
}

class CachingTokenEventLogger(cacheEntryTTL: FiniteDuration) extends TokenEventLogger {

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


  private val backendCache = CacheBuilder.newBuilder()
    .expireAfterWrite(cacheEntryTTL._1, cacheEntryTTL._2)
    .maximumSize(10000)
    .build[String, Object]()

  override def outOfTokens(backend: String): Unit = {
    backendCache.put(backend, new Object())
  }

  override def getLimitedBackends: Set[String] = {
    import scala.collection.JavaConverters._
    backendCache.asMap().keySet().asScala.toSet
  }

}
