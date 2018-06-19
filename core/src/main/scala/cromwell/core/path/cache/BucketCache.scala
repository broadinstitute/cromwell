package cromwell.core.path.cache

import com.google.common.cache.Cache

/**
  * Use to cache meta information about buckets. Keys are bucket names.
  */
abstract class BucketCache[A <: Object](cache: Cache[String, A]) extends CacheHelper(cache)
