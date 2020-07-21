package cromwell.backend.standard.callcaching
import cats.data.NonEmptyList
import cromwell.backend.BackendCacheHitCopyingActor.CopyOutputsCommand
import cromwell.core.io.{IoCommand, IoCopyCommand}
import cromwell.services.CallCaching.CallCachingEntryId


object CopyingActorBlacklistCacheSupport {
  trait HasFormatting {
    def metricFormat: String = getClass.getName.toLowerCase.split('$').last
  }

  sealed trait Verb extends HasFormatting
  case object Read extends Verb
  case object Write extends Verb

  sealed trait EntityType extends HasFormatting
  case object Hit extends EntityType
  case object Bucket extends EntityType

  sealed trait CacheReadType
  case object ReadHitOnly
  case object ReadHitAndBucket
}

trait CopyingActorBlacklistCacheSupport {
  this: StandardCacheHitCopyingActor =>

  import CopyingActorBlacklistCacheSupport._

  def handleBlacklistingForGenericFailure(): Unit = {
    // Not a forbidden failure so do not blacklist the bucket but do blacklist the hit.
    for {
      data <- stateData
      cache <- standardParams.blacklistCache
      _ = blacklistAndMetricHit(cache, data.cacheHit)
    } yield ()
    ()
  }

  /* Whitelist by bucket and hit if appropriate. */
  def handleWhitelistingForSuccess(command: IoCommand[_]): Unit = {
    for {
      cache <- standardParams.blacklistCache
      data <- stateData
      _ = whitelistAndMetricHit(cache, data.cacheHit)
      copy <- Option(command) collect { case c: IoCopyCommand => c }
      prefix <- extractBlacklistPrefix(copy.source.toString)
      _ = whitelistAndMetricBucket(cache, prefix)
    } yield ()
    ()
  }

  def publishBlacklistMetric(blacklistCache: BlacklistCache, verb: Verb, entityType: EntityType, value: BlacklistStatus): Unit = {
    val group = blacklistCache.name.getOrElse("none")
    val metricPath = NonEmptyList.of(
      "job",
      "callcaching", "blacklist", verb.metricFormat, entityType.metricFormat, group, value.toString)
    increment(metricPath)
  }

  def blacklistAndMetricHit(blacklistCache: BlacklistCache, hit: CallCachingEntryId): Unit = {
    blacklistCache.getBlacklistStatus(hit) match {
      case UntestedCacheResult =>
        blacklistCache.blacklist(hit)
        publishBlacklistMetric(blacklistCache, Write, Hit, value = BadCacheResult)
      case BadCacheResult =>
      // Not a surprise, race conditions abound in cache hit copying. Do not overwrite with the same value or
      // multiply publish metrics for this hit.
      case GoodCacheResult =>
        // This hit was thought to be good but now a copy has failed for permissions reasons. Be conservative and
        // mark the hit as BadCacheResult and log this strangeness.
        log.warning(
          "Cache hit {} found in GoodCacheResult blacklist state, but cache hit copying has failed for permissions reasons. Overwriting status to BadCacheResult state.",
          hit.id)
        blacklistCache.blacklist(hit)
        publishBlacklistMetric(blacklistCache, Write, Hit, value = BadCacheResult)
    }
  }

  def blacklistAndMetricBucket(blacklistCache: BlacklistCache, bucket: String): Unit = {
    blacklistCache.getBlacklistStatus(bucket) match {
      case UntestedCacheResult =>
        blacklistCache.blacklist(bucket)
        publishBlacklistMetric(blacklistCache, Write, Bucket, value = BadCacheResult)
      case BadCacheResult =>
      // Not a surprise, race conditions abound in cache hit copying. Do not overwrite with the same value or
      // multiply publish metrics for this bucket.
      case GoodCacheResult =>
        // This bucket was thought to be good but now a copy has failed for permissions reasons. Be conservative and
        // mark the bucket as BadCacheResult and log this strangeness.
        log.warning(
          "Bucket {} found in GoodCacheResult blacklist state, but cache hit copying has failed for permissions reasons. Overwriting status to BadCacheResult state.",
          bucket)
        blacklistCache.blacklist(bucket)
        publishBlacklistMetric(blacklistCache, Write, Bucket, value = BadCacheResult)
    }
  }

  def whitelistAndMetricHit(blacklistCache: BlacklistCache, hit: CallCachingEntryId): Unit = {
    blacklistCache.getBlacklistStatus(hit) match {
      case UntestedCacheResult =>
        blacklistCache.whitelist(hit)
        publishBlacklistMetric(blacklistCache, Write, Hit, value = GoodCacheResult)
      case GoodCacheResult => // This hit is already known to be good, no need to rewrite or spam metrics.
      case BadCacheResult =>
        // This is surprising, a hit that we failed to copy before has now been the source of a successful copy.
        // Don't overwrite this to GoodCacheResult, hopefully there are less weird cache hits out there.
        log.warning(
          "Cache hit {} found in BadCacheResult blacklist state, not overwriting to GoodCacheResult despite successful copy.",
          hit.id)
    }
  }

  def whitelistAndMetricBucket(blacklistCache: BlacklistCache, bucket: String): Unit = {
    blacklistCache.getBlacklistStatus(bucket) match {
      case UntestedCacheResult =>
        blacklistCache.whitelist(bucket)
        publishBlacklistMetric(blacklistCache, Write, Bucket, value = GoodCacheResult)
      case GoodCacheResult => // This bucket is already known to be good, no need to rewrite or spam metrics.
      case BadCacheResult =>
        // This is surprising, a bucket that we failed to copy from before for auth reasons has now been the source
        // of a successful copy. Don't overwrite this to GoodCacheResult, hopefully there are less weird cache hits out there.
        log.warning(
          "Bucket {} found in BadCacheResult blacklist state, not overwriting to GoodCacheResult despite successful copy.",
          bucket)
    }
  }

  def publishBlacklistReadMetrics(command: CopyOutputsCommand, cacheHit: CallCachingEntryId, cacheReadType: Product) = {
    for {
      c <- standardParams.blacklistCache
      hitBlacklistStatus = c.getBlacklistStatus(cacheHit)
      // If blacklisting is on the hit cache is always checked so publish a hit read metric.
      _ = publishBlacklistMetric(c, Read, Hit, hitBlacklistStatus)
      // Conditionally publish the bucket read if the backend supports bucket / prefix blacklisting and the bucket was read.
      _ <- Option(cacheReadType).collect { case ReadHitAndBucket => () }
      path = sourcePathFromCopyOutputsCommand(command)
      prefix <- extractBlacklistPrefix(path)
      bucketBlacklistStatus = c.getBlacklistStatus(prefix)
      _ = publishBlacklistMetric(c, Read, Bucket, bucketBlacklistStatus)
    } yield ()
  }

  def isSourceBlacklisted(command: CopyOutputsCommand): Boolean = {
    val path = sourcePathFromCopyOutputsCommand(command)
    (for {
      cache <- standardParams.blacklistCache
      prefix <- extractBlacklistPrefix(path)
      value = cache.getBlacklistStatus(prefix)
    } yield value == BadCacheResult).getOrElse(false)
  }

  def isSourceBlacklisted(hit: CallCachingEntryId): Boolean = {
    (for {
      cache <- standardParams.blacklistCache
      value = cache.getBlacklistStatus(hit)
    } yield value == BadCacheResult).getOrElse(false)
  }
}
