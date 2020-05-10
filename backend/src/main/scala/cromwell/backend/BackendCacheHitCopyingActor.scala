package cromwell.backend

import cromwell.backend.MetricableCacheCopyErrorCategory.MetricableCacheCopyErrorCategory
import cromwell.core.JobKey
import cromwell.core.simpleton.WomValueSimpleton
import cromwell.services.CallCaching.CallCachingEntryId

object BackendCacheHitCopyingActor {
  final case class CopyOutputsCommand(womValueSimpletons: Seq[WomValueSimpleton], jobDetritusFiles: Map[String, String], cacheHit: CallCachingEntryId, returnCode: Option[Int])

  // 'declined' means that the copying actor declined to attempt copying the cache hit. Currently this is true only
  // for blacklisting reasons.
  final case class CopyingOutputsFailedResponse(jobKey: JobKey, cacheCopyAttempt: Int, failure: CacheCopyError, declined: Boolean = false)

  sealed trait CacheCopyError
  final case class LoggableCacheCopyError(failure: Throwable) extends CacheCopyError
  final case class MetricableCacheCopyError(failureCategory: MetricableCacheCopyErrorCategory) extends CacheCopyError
}

object MetricableCacheCopyErrorCategory {
  sealed trait MetricableCacheCopyErrorCategory {
    // The 'stripSuffix is necessary to handle scala 'case object's which get '$' added to the class name during java translation.
    override def toString: String = getClass.getSimpleName.stripSuffix("$").toLowerCase
  }
  final case object BucketBlacklisted extends MetricableCacheCopyErrorCategory
  final case object HitBlacklisted extends MetricableCacheCopyErrorCategory
}
