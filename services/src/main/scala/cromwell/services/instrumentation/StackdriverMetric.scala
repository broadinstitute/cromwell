package cromwell.services.instrumentation

import scala.concurrent.duration.FiniteDuration

sealed trait StackdriverMetric

case class StackdriverTiming(bucket: CromwellBucket, value: FiniteDuration) extends StackdriverMetric
