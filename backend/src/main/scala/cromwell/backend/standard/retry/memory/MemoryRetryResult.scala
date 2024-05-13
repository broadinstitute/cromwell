package cromwell.backend.standard.retry.memory

import common.validation.Validation.MemoryRetryMultiplierRefined

/**
  * Result of evaluating an attempt as a memory-retry candidate, encapsulating instructions
  * for configuring the next attempt if applicable.
  *
  * @param oomDetected        Did the previous attempt OOM?
  * @param factor             User-configured factor
  * @param previousMultiplier Multiplier used for the previous attempt
  */
case class MemoryRetryResult(oomDetected: Boolean,
                             factor: Option[MemoryRetryMultiplierRefined],
                             previousMultiplier: Option[Double]
)

case object MemoryRetryResult {
  def none: MemoryRetryResult = MemoryRetryResult(oomDetected = false, None, None)
}
