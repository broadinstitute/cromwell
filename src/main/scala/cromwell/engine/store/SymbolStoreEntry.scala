package cromwell.engine.store

import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.WdlValue
import cromwell.engine.store.SymbolStore.SymbolStoreKey

object SymbolStoreEntry {
  private def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
    val lastIndex = fullyQualifiedName.lastIndexOf(".")
    (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
  }

  def apply(fullyQualifiedName: FullyQualifiedName, wdlValue: WdlValue, input: Boolean): SymbolStoreEntry = {
    val (scope, name) = splitFqn(fullyQualifiedName)
    val key = SymbolStoreKey(scope, name, iteration = None, input = true)
    SymbolStoreEntry(key, wdlValue.wdlType, Some(wdlValue))
  }
}

case class SymbolStoreEntry(key: SymbolStoreKey, wdlType: WdlType, wdlValue: Option[WdlValue]) {
  def isInput: Boolean = key.input
  def isOutput: Boolean = !isInput
  def scope: String = key.scope
}