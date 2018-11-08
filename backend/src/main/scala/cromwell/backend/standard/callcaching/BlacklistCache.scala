package cromwell.backend.standard.callcaching

import com.google.common.cache.{CacheBuilder, CacheLoader}
import cromwell.core.CacheConfig

case class BlacklistCache(config: CacheConfig) {
  val cache = {
    // Queries to the blacklist cache return false by default (i.e. not blacklisted).
    val falseLoader = new CacheLoader[String, java.lang.Boolean]() {
      override def load(key: String): java.lang.Boolean = false
    }

    CacheBuilder.
      newBuilder().
      concurrencyLevel(config.concurrency).
      maximumSize(config.size).
      expireAfterWrite(config.ttl.length, config.ttl.unit).
      build[String, java.lang.Boolean](falseLoader)
  }

  def isBlacklisted(bucket: String): Boolean = cache.get(bucket)

  def blacklist(bucket: String): Unit = cache.put(bucket, true)
}
