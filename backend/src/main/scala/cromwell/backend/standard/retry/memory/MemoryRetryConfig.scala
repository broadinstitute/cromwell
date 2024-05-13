package cromwell.backend.standard.retry.memory

import common.validation.Validation.MemoryRetryMultiplierRefined

/**
  * Result of evaluating an attempt as a memory-retry candidate, encapsulating instructions
  * for configuring the next attempt if applicable.
  *
  * @param oomDetected        Did the last attempt OOM?
  * @param factor             User-configured factor
  * @param previousMultiplier Multiplier used for the last attempt
  */
case class MemoryRetryConfig(oomDetected: Boolean,
                             factor: Option[MemoryRetryMultiplierRefined],
                             previousMultiplier: Option[Double])

case object MemoryRetryConfig {
  def none: MemoryRetryConfig = MemoryRetryConfig(false, None, None)
}
