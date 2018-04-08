package cromwell.backend.impl.jes

import cromwell.services.keyvalue.KeyValueServiceActor._
import common.validation.ErrorOr.ErrorOr
import cats.syntax.validated._
import cats.syntax.apply._

import scala.util.{Failure, Success, Try}
import JesBackendLifecycleActorFactory.preemptionCountKey
import JesBackendLifecycleActorFactory.failureCountKey

case class PreviousRetryReasons(preempted: Int, failed: Int)

object PreviousRetryReasons {
  def tryApply(prefetchedKvEntries: Map[String, KvResponse], attemptNumber: Int): ErrorOr[PreviousRetryReasons] = {
    val validatedPreemptionCount = validatedKvResponse(prefetchedKvEntries.get(preemptionCountKey), preemptionCountKey)
    val validatedFailedCount = validatedKvResponse(prefetchedKvEntries.get(failureCountKey), failureCountKey)

    (validatedPreemptionCount, validatedFailedCount) mapN { PreviousRetryReasons.apply }
  }

  def apply(knownPreemptedCount: Int, knownFailureCount: Int, attempt: Int): PreviousRetryReasons = {
    // If we have anything unaccounted for, we can top up the general retry count.
    // NB: 'attempt' is 1-indexed, so, magic number:
    // NB2: for sanity's sake, I won't let this unaccounted for drop below 0, just in case...
    val unaccountedFor = Math.max(attempt - 1 - knownPreemptedCount - knownFailureCount, 0)
    PreviousRetryReasons(knownPreemptedCount, knownFailureCount + unaccountedFor)
  }

  private def validatedKvResponse(r: Option[KvResponse], fromKey: String): ErrorOr[Int] = r match {
    case Some(KvPair(_, v)) => validatedInt(v, fromKey)
    case Some(_: KvKeyLookupFailed) => 0.validNel
    case Some(KvFailure(_, failure)) => s"Failed to get key $fromKey: ${failure.getMessage}".invalidNel
    case Some(_: KvPutSuccess) => s"Programmer Error: Got a KvPutSuccess from a Get request...".invalidNel
    case None => s"Programmer Error: Engine made no effort to prefetch $fromKey".invalidNel
  }

  private def validatedInt(s: String, fromKey: String): ErrorOr[Int] = {
    Try(s.toInt) match {
      case Success(i) => i.validNel
      case Failure(_) => s"Unexpected value found in the KV store: $fromKey='$s'".invalidNel
    }
  }
}
