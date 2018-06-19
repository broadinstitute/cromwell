package cromwell.filesystems.gcs.cache

import cats.effect.IO
import com.google.cloud.storage.Storage
import com.google.common.cache.{Cache, CacheBuilder}
import cromwell.filesystems.gcs.cache.GcsBucketCache.GcsBucketInformation

/**
  * Simple trait with a function returning information about a bucket, possibly asynchronously.
  */
trait GcsBucketInformationStore {
  def information(bucket: String): IO[GcsBucketInformation]
}

/**
  * Defines behavior on how meta information about buckets can be retrieved and stored.
  */
object GcsBucketInformationPolicies {
  sealed trait GcsBucketInformationPolicy {
    /**
      * Create a store corresponding to the current policy
      */
    def toStore(cloudStorage: Storage, projectId: String): GcsBucketInformationStore
  }

  /**
    * Use to return a default static bucket information value
    */
  case object Default extends GcsBucketInformationPolicy {
    private val value = GcsBucketInformation(requesterPays = false)
    private val store = new GcsBucketInformationStore {
      def information(bucket: String): IO[GcsBucketInformation] = IO.pure(value)
    }
    override def toStore(cloudStorage: Storage, projectId: String) = store
  }

  /**
    * Use to force retrieval of new bucket information for every path being built.
    */
  case object NoCache extends GcsBucketInformationPolicy {
    private val noCache = CacheBuilder.newBuilder().maximumSize(0).build[String, GcsBucketInformation]()
    override def toStore(cloudStorage: Storage, projectId: String) = new GcsBucketCache(cloudStorage, noCache, projectId) with GcsBucketInformationStore {
      override def information(bucket: String) = getCachedValue(bucket)
    }
  }

  /**
    * Use to cache bucket information in the provided Guava cache
    */
  case class CacheInformation(cache: Cache[String, GcsBucketInformation]) extends GcsBucketInformationPolicy {
    override def toStore(cloudStorage: Storage, projectId: String) = new GcsBucketCache(cloudStorage, cache, projectId) with GcsBucketInformationStore {
      override def information(bucket: String) = getCachedValue(bucket)
    }
  }
}
