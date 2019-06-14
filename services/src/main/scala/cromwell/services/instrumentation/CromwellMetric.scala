package cromwell.services.instrumentation

import cats.data.NonEmptyList
import cromwell.services.instrumentation.CromwellInstrumentation.InstrumentationPath

import scala.concurrent.duration.FiniteDuration

object CromwellBucket {
  def apply(singleString: String): CromwellBucket = CromwellBucket(List.empty, NonEmptyList.of(singleString))
}

/**
  * Represents an instrumentation bucket. Separate prefix from path.
  * For example a bucket instrumenting number of workflow submitted would would have a "workflow" prefix and an  "submitted" key
  * This allow for deferred insertion of control elements in between prefix and key when building a bucket string.
  */
case class CromwellBucket(prefix: List[String], path: InstrumentationPath) {
  def expand(key: String): CromwellBucket = this.copy(path = path.concatNel(NonEmptyList.of(key)))
}

/**
  * Generic type for Cromwell Instrumentation Metrics
  */
sealed trait CromwellMetric

object CromwellCount {
  def apply(bucket: CromwellBucket, value: Long): CromwellCount = new CromwellCount(bucket, value)
  def unapply(count: CromwellCount): Option[(CromwellBucket, Long, Double)] = Option((count.bucket, count.value, count.sampling))
}
/**
  * Count occurrences of an event
  */
class CromwellCount(val bucket: CromwellBucket, val value: Long, val sampling: Double = 1.0D) extends CromwellMetric {
  override def toString: String = s"CromwellCount($bucket, $value, $sampling)"
}

/**
  * Convenience class to increment a count (+1)
  */
case class CromwellIncrement(override val bucket: CromwellBucket) extends CromwellCount(bucket, 1)

/**
  * Measures a time value
  */
case class CromwellTiming(bucket: CromwellBucket, value: FiniteDuration, sampling: Double = 1.0D) extends CromwellMetric

/**
  * Measures a gauge value
  */
case class CromwellGauge(bucket: CromwellBucket, value: Long) extends CromwellMetric
