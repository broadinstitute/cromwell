package cromwell.engine

import cromwell.binding.{Call, WdlBinding}

/**
 * Corresponds to the "execution table" of our discussions.
 */
class ExecutionStatusStore(binding: WdlBinding) {
  /** FIXME this implementation only sorta works for "hello world" and is
    * FIXME not adequate for Sprint 2/3. */
  def runnableCalls: Set[Call] = binding.workflow.calls
}
