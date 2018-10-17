package cromwell.backend.impl.tes

/**
  * GRANULAR: The TES payload will list outputs individually in the outputs section of the payload
  * ROOT: The TES payload will list only the root execution directory in the outputs section of the payload
  */
object OutputMode extends Enumeration {
  type OutputMode = Value
  val GRANULAR, ROOT = Value
}
