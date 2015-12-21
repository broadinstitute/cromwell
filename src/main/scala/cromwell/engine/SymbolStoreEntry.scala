package cromwell.engine

import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.WdlValue

object SymbolStoreEntry {
  private def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
    val lastIndex = fullyQualifiedName.lastIndexOf(".")
    (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
  }

  def apply(fullyQualifiedName: FullyQualifiedName, wdlValue: WdlValue, symbolHash: SymbolHash, input: Boolean): SymbolStoreEntry = {
    val (scope, name) = splitFqn(fullyQualifiedName)
    val key = SymbolStoreKey(scope, name, index = None, input)
    SymbolStoreEntry(key, wdlValue.wdlType, Option(wdlValue), Option(symbolHash))
  }

  def toWorkflowOutputs(t: Traversable[SymbolStoreEntry]): WorkflowOutputs = t.map { e =>
    s"${e.key.scope}.${e.key.name}" -> CallOutput(e.wdlValue.get, e.symbolHash.get)
  }.toMap

  def toCallOutputs(traversable: Traversable[SymbolStoreEntry]): CallOutputs = traversable.map { entry =>
    entry.key.name -> CallOutput(entry.wdlValue.get, entry.symbolHash.get)
  }.toMap
}

case class SymbolStoreEntry(key: SymbolStoreKey, wdlType: WdlType, wdlValue: Option[WdlValue], symbolHash: Option[SymbolHash]) {
  def isInput: Boolean = key.input
  def isOutput: Boolean = !isInput
  def scope: String = key.scope
}

case class SymbolStoreKey(scope: String, name: String, index: Option[Int], input: Boolean) {
  def fqn: FullyQualifiedName = s"$scope.$name"
}