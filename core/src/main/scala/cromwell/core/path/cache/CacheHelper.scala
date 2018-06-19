package cromwell.core.path.cache

import cats.effect.IO
import com.google.common.cache.Cache

/**
  * Generic class providing methods to help with caching.
  */
abstract class CacheHelper[K <: Object, V <: Object](cache: Cache[K, V]) {

  /**
    * Given a key K, returns value V associated with it.
    * This will only be called when the key is not in the cache.
    * This method should ONLY be called from this trait.
    */
  protected def retrieve(key: K): IO[V]

  private def getOrInsert(key: K) = Option(cache.getIfPresent(key))
    .map(IO.pure)
    .getOrElse({
      retrieve(key) map { value =>
        cache.put(key, value)
        value
      }
    })

  /**
    * Returns the value associated with the key.
    * Use the cache if possible, otherwise call [[retrieve]]
    */
  final def getCachedValue(key: K): IO[V] = getOrInsert(key)
}
