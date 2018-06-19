package cromwell.core.path.cache

import java.nio.file.FileSystem

import com.google.common.cache.Cache

/**
  * Use to cache filesystems by associated bucket.
  * Be careful to use different caches for different credentials. 
  */
abstract class FileSystemCache[A <: FileSystem](cache: Cache[String, A]) extends CacheHelper(cache)
