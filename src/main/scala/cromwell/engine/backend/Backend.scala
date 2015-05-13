package cromwell.engine.backend

import cromwell.binding.WdlExpression.ScopedLookupFunction
import cromwell.binding.values.WdlValue
import cromwell.binding.{Call, TaskOutput}
import cromwell.engine.SymbolStore

import scala.util.Try

/**
 * trait to be implemented by concrete backends.
 */
trait Backend {

  /** Curried version that a caller is expected to fix for a Call and symbol store. */
  private def lookupFunction(call: Call, symbolStore: SymbolStore)(localName: String): WdlValue =
    symbolStore.locallyQualifiedInputs(call).find { _._1 == localName }.get._2

  def scopedLookupFunction(call: Call, symbolStore: SymbolStore): ScopedLookupFunction =
    lookupFunction(call, symbolStore)

  /**
   * Execute the specified command line using the provided symbol store, evaluating the task outputs to produce
   * a mapping of local task output names to WDL values.
   */
  def executeCommand(commandLine: String, call: Call, taskOutputs: Set[TaskOutput], symbolStore: SymbolStore): Map[String, Try[WdlValue]]

}
