package cromwell.backend.standard.callcaching

import akka.event.LoggingAdapter
import com.google.common.cache.{CacheBuilder, CacheLoader, LoadingCache}
import com.typesafe.config.Config
import cromwell.core.{CacheConfig, HasWorkflowIdAndSources}
import mouse.boolean._
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps

object CallCachingBlacklistManager {
  object Defaults {
    object Groupings {
      val Concurrency = 10000
      val Size = 1000L
      val Ttl = 2 hours
    }
    object Hits {
      val Concurrency = 10000
      val Size = 20000L
      val Ttl = 1 hour
    }
    object Buckets {
      val Concurrency = 10000
      val Size = 1000L
      val Ttl = 1 hour
    }
  }
}

class CallCachingBlacklistManager(rootConfig: Config, logger: LoggingAdapter) {

  // Defined if "call-caching.blacklist-cache.enabled = true".
  private val blacklistCacheConfig: Option[Unit] =
    rootConfig.getOrElse("call-caching.blacklist-cache.enabled", false).option(())

  // Defined if `blacklistCacheConfig` is defined and "call-caching.blacklist-cache.groupings.workflow-option" is defined.
  private val blacklistGroupingWorkflowOptionKey: Option[String] = for {
    _ <- blacklistCacheConfig // Only return a groupings cache if blacklisting is enabled.
    workflowOption <- rootConfig.as[Option[String]]("call-caching.blacklist-cache.groupings.workflow-option")
  } yield workflowOption

  // Defined if `blacklistGroupingWorkflowOptionKey` is defined.
  private val blacklistGroupingCacheConfig: Option[CacheConfig] = {
    import CallCachingBlacklistManager.Defaults.Groupings._
    for {
      _ <- blacklistGroupingWorkflowOptionKey
      groupingsOption = rootConfig.as[Option[Config]] ("call-caching.blacklist-cache.groupings")
      conf = CacheConfig.config(groupingsOption, defaultConcurrency = Concurrency, defaultSize = Size, defaultTtl = Ttl)
    } yield conf
  }

  // Defined if `blacklistCacheConfig` is defined.
  private val blacklistBucketCacheConfig: Option[CacheConfig] = {
    import CallCachingBlacklistManager.Defaults.Buckets._
    for {
      _ <- blacklistCacheConfig
      bucketsOption = rootConfig.as[Option[Config]]("call-caching.blacklist-cache.buckets")
      conf = CacheConfig.config(bucketsOption, defaultConcurrency = Concurrency, defaultSize = Size, defaultTtl = Ttl)
    } yield conf
  }

  // Defined if `blacklistCacheConfig` is defined.
  private val blacklistHitCacheConfig: Option[CacheConfig] = {
    import CallCachingBlacklistManager.Defaults.Hits._
    for {
      _ <- blacklistCacheConfig
      hitsOption = rootConfig.as[Option[Config]]("call-caching.blacklist-cache.hits")
      conf = CacheConfig.config(hitsOption, defaultConcurrency = Concurrency, defaultSize = Size, defaultTtl = Ttl)
    } yield conf
  }

  // If configuration allows, build a cache of blacklist groupings to BlacklistCaches.
  private val blacklistGroupingsCache: Option[LoadingCache[String, BlacklistCache]] = {
    def buildBlacklistGroupingsCache(groupingConfig: CacheConfig, bucketConfig: CacheConfig, hitConfig: CacheConfig): LoadingCache[String, BlacklistCache] = {
      val emptyBlacklistCacheLoader = new CacheLoader[String, BlacklistCache]() {
        override def load(key: String): BlacklistCache = new GroupingBlacklistCache(
          bucketCacheConfig = bucketConfig,
          hitCacheConfig = hitConfig,
          group = key
        )
      }

      CacheBuilder.
        newBuilder().
        concurrencyLevel(groupingConfig.concurrency).
        maximumSize(groupingConfig.size).
        expireAfterWrite(groupingConfig.ttl.length, groupingConfig.ttl.unit).
        build[String, BlacklistCache](emptyBlacklistCacheLoader)
    }

    for {
      groupingsConfig <- blacklistGroupingCacheConfig
      bucketsConfig <- blacklistBucketCacheConfig
      hitsConfig <- blacklistHitCacheConfig
    } yield buildBlacklistGroupingsCache(groupingsConfig, bucketsConfig, hitsConfig)
  }

  /**
    * If configured return a group blacklist cache, otherwise if configured return a root workflow cache,
    * otherwise return nothing.
    */
  def blacklistCacheFor(workflow: HasWorkflowIdAndSources): Option[BlacklistCache] = {
    // If configuration is set up for blacklist groups and a blacklist group is specified in workflow options,
    // get the BlacklistCache for the group.
    val groupBlacklistCache: Option[BlacklistCache] = for {
      groupings <- blacklistGroupingsCache
      groupKey <- blacklistGroupingWorkflowOptionKey
      groupFromWorkflowOptions <- workflow.sources.workflowOptions.get(groupKey).toOption
    } yield groupings.get(groupFromWorkflowOptions)

    // Build a blacklist cache for a single, ungrouped root workflow.
    def rootWorkflowBlacklistCache: Option[BlacklistCache] = for {
      bucketConfig <- blacklistBucketCacheConfig
      hitConfig <- blacklistHitCacheConfig
    } yield new RootWorkflowBlacklistCache(bucketCacheConfig = bucketConfig, hitCacheConfig = hitConfig)

    // Return the group blacklist cache if available, otherwise a blacklist cache for the root workflow.
    val maybeCache = groupBlacklistCache orElse rootWorkflowBlacklistCache
    maybeCache collect {
      case group: GroupingBlacklistCache =>
        logger.info("Workflow {} using group blacklist cache '{}' containing blacklist status for {} hits and {} buckets.",
          workflow.id, group.group, group.hitCache.size(), group.bucketCache.size())
      case _: RootWorkflowBlacklistCache =>
        logger.info("Workflow {} using root workflow blacklist cache.", workflow.id)
    }
    maybeCache
  }
}
