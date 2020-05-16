package cromwell.engine.workflow

import akka.event.LoggingAdapter
import com.google.common.cache.{CacheBuilder, CacheLoader, LoadingCache}
import com.typesafe.config.Config
import cromwell.backend.standard.callcaching.{BlacklistCache, GroupBlacklistCache, RootWorkflowBlacklistCache}
import cromwell.core.CacheConfig
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps

class CallCachingBlacklistManager(rootConfig: Config, logger: LoggingAdapter) {

  // Defined only if "call-caching.blacklist-cache" is defined in config and "enabled = true".
  private val blacklistCacheConfig: Option[CacheConfig] = {
    for {
      conf <- rootConfig.as[Option[Config]] ("call-caching.blacklist-cache")
      cacheConfig <- CacheConfig.optionalConfig(conf, defaultConcurrency = 1000, defaultSize = 1000, defaultTtl = 2 hours)
    } yield cacheConfig
  }

  // Defined only if `blacklistCacheConfig` is defined and "call-caching.blacklist-cache.groupings.workflow-option" is defined.
  private val blacklistGroupingWorkflowOptionKey: Option[String] = for {
    _ <- blacklistCacheConfig // Only return a groupings cache if blacklisting is enabled.
    workflowOption <- rootConfig.as[Option[String]]("call-caching.blacklist-cache.groupings.workflow-option")
  } yield workflowOption

  // Defined only if `blacklistGroupingWorkflowOptionKey` is defined.
  private val blacklistGroupingCacheConfig: Option[CacheConfig] = {
    for {
      _ <- blacklistGroupingWorkflowOptionKey
      groupingsOption = rootConfig.as[Option[Config]] ("call-caching.blacklist-cache.groupings")
      conf = CacheConfig.config(groupingsOption, defaultConcurrency = 1000, defaultSize = 1000, defaultTtl = 2 hours)
    } yield conf
  }

  // Defined if `blacklistCacheConfig` is defined.
  private val blacklistBucketCacheConfig: Option[CacheConfig] = for {
    _ <- blacklistCacheConfig
    bucketsOption = rootConfig.as[Option[Config]]("call-caching.blacklist-cache.buckets")
    conf = CacheConfig.config(bucketsOption, defaultConcurrency = 1000, defaultSize = 1000, defaultTtl = 1 hour)
  } yield conf

  // Defined if `blacklistCacheConfig` is defined.
  private val blacklistHitCacheConfig: Option[CacheConfig] = for {
    _ <- blacklistCacheConfig
    hitsOption = rootConfig.as[Option[Config]]("call-caching.blacklist-cache.hits")
    conf = CacheConfig.config(hitsOption, defaultConcurrency = 1000, defaultSize = 50000, defaultTtl = 1 hour)
  } yield conf

  // If configuration allows, build a cache of blacklist groupings to BlacklistCaches.
  private val blacklistGroupingsCache: Option[LoadingCache[String, BlacklistCache]] = {
    def buildBlacklistGroupingsCache(groupingConfig: CacheConfig, bucketConfig: CacheConfig, hitConfig: CacheConfig): LoadingCache[String, BlacklistCache] = {
      val emptyBlacklistCacheLoader = new CacheLoader[String, BlacklistCache]() {
        override def load(key: String): BlacklistCache = new GroupBlacklistCache(
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
  def blacklistCacheFor(workflow: workflowstore.WorkflowToStart): Option[BlacklistCache] = {
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
      case group: GroupBlacklistCache =>
        logger.info("Workflow {} using group blacklist cache '{}' containing blacklist status for {} hits and {} buckets.",
          workflow.id, group.group, group.hitCache.size(), group.bucketCache.size())
      case _: RootWorkflowBlacklistCache =>
        logger.info("Workflow {} using root workflow blacklist cache.", workflow.id)
    }
    maybeCache
  }
}
