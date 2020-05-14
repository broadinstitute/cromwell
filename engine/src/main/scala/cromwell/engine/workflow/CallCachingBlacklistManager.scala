package cromwell.engine.workflow

import com.google.common.cache.{CacheBuilder, CacheLoader, LoadingCache}
import com.typesafe.config.Config
import cromwell.backend.standard.callcaching.{BlacklistCache, GroupedBlacklistCache, RootWorkflowBlacklistCache}
import cromwell.core.CacheConfig
import net.ceedubs.ficus.Ficus._

import scala.concurrent.duration._
import scala.language.postfixOps

class CallCachingBlacklistManager(config: Config) {
  private val blacklistCacheConfig: Option[Config] = config.as[Option[Config]]("call-caching.blacklist-cache")

  private val blacklistGroupCacheConfig: Option[Config] = blacklistCacheConfig.flatMap(_.as[Option[Config]]("groupings"))

  private val blacklistGroupWorkflowOptionKey: Option[String] = blacklistGroupCacheConfig.flatMap(_.as[Option[String]]("workflow-option"))

  // If configuration allows, build a cache of blacklist groupings to BlacklistCaches.
  private val blacklistGroupingsCache: Option[LoadingCache[String, BlacklistCache]] = {
    def buildBlacklistGroupingsCache(cacheConfig: CacheConfig): LoadingCache[String, BlacklistCache] = {
      val emptyBlacklistCacheLoader = new CacheLoader[String, BlacklistCache]() {
        override def load(key: String): BlacklistCache = new GroupedBlacklistCache(cacheConfig, group = key)
      }

      CacheBuilder.
        newBuilder().
        concurrencyLevel(cacheConfig.concurrency).
        maximumSize(cacheConfig.size).
        expireAfterWrite(cacheConfig.ttl.length, cacheConfig.ttl.unit).
        build[String, BlacklistCache](emptyBlacklistCacheLoader)
    }

    for {
      conf <- blacklistGroupCacheConfig
      _ <- blacklistGroupWorkflowOptionKey // Only build groupings if a workflow option key is specified in config.
      groupingsConfig <- CacheConfig.optionalConfig(conf, defaultConcurrency = 10000, defaultSize = 1000, defaultTtl = 2 hours)
    } yield buildBlacklistGroupingsCache(groupingsConfig)
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
      groupKey <- blacklistGroupWorkflowOptionKey
      groupFromWorkflowOptions <- workflow.sources.workflowOptions.get(groupKey).toOption
    } yield groupings.get(groupFromWorkflowOptions)

    // Build a blacklist cache for a single, ungrouped root workflow.
    def rootWorkflowBlacklistCache: Option[BlacklistCache] = {
      for {
        config <- blacklistCacheConfig
        cacheConfig <- CacheConfig.optionalConfig(config, defaultConcurrency = 1000, defaultSize = 1000, defaultTtl = 1 hour)
      } yield new RootWorkflowBlacklistCache(cacheConfig)
    }

    // Return the group blacklist cache if available, otherwise a blacklist cache for the root workflow.
    groupBlacklistCache orElse rootWorkflowBlacklistCache
  }
}
