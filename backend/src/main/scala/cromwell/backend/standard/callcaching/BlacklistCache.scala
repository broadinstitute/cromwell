package cromwell.backend.standard.callcaching

import com.google.common.cache.{CacheBuilder, CacheLoader}
import cromwell.core.CacheConfig
import cromwell.services.CallCaching.CallCachingEntryId

sealed trait BlacklistStatus
case object KnownBad extends BlacklistStatus
case object KnownGood extends BlacklistStatus
case object Unknown extends BlacklistStatus

sealed abstract class BlacklistCache(val bucketCacheConfig: CacheConfig,
                                     val hitCacheConfig: CacheConfig,
                                     val name: Option[String]) {
  val bucketCache = {
    // Queries to the bucket blacklist cache return Unknown by default.
    val unknownLoader = new CacheLoader[String, BlacklistStatus]() {
      override def load(key: String): BlacklistStatus = Unknown
    }

    CacheBuilder.
      newBuilder().
      concurrencyLevel(bucketCacheConfig.concurrency).
      maximumSize(bucketCacheConfig.size).
      expireAfterWrite(bucketCacheConfig.ttl.length, bucketCacheConfig.ttl.unit).
      build[String, BlacklistStatus](unknownLoader)
  }

  val hitCache = {
    // Queries to the hit blacklist cache return false by default (i.e. not blacklisted).
    val unknownLoader = new CacheLoader[CallCachingEntryId, BlacklistStatus]() {
      override def load(key: CallCachingEntryId): BlacklistStatus = Unknown
    }

    CacheBuilder.
      newBuilder().
      concurrencyLevel(hitCacheConfig.concurrency).
      maximumSize(hitCacheConfig.size).
      expireAfterWrite(hitCacheConfig.ttl.length, hitCacheConfig.ttl.unit).
      build[CallCachingEntryId, BlacklistStatus](unknownLoader)
  }

  def isBlacklisted(hit: CallCachingEntryId): BlacklistStatus = hitCache.get(hit)

  def isBlacklisted(bucket: String): BlacklistStatus = bucketCache.get(bucket)

  def blacklistHit(hit: CallCachingEntryId): Unit = hitCache.put(hit, KnownBad)

  def blacklistBucket(bucket: String): Unit = bucketCache.put(bucket, KnownBad)

  def whitelistHit(hit: CallCachingEntryId): Unit = hitCache.put(hit, KnownGood)

  def whitelistBucket(bucket: String): Unit = bucketCache.put(bucket, KnownGood)
}

class RootWorkflowBlacklistCache(bucketCacheConfig: CacheConfig, hitCacheConfig: CacheConfig) extends
  BlacklistCache(bucketCacheConfig = bucketCacheConfig, hitCacheConfig = hitCacheConfig, name = None)

class GroupBlacklistCache(bucketCacheConfig: CacheConfig, hitCacheConfig: CacheConfig, val group: String) extends
  BlacklistCache(bucketCacheConfig = bucketCacheConfig, hitCacheConfig = hitCacheConfig, name = Option(group))
