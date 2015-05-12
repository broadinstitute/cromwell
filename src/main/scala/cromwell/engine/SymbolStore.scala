package cromwell.engine

import cromwell.binding.types.WdlType
import cromwell.binding._
import cromwell.binding.values.WdlValue

import scala.collection.mutable


case class SymbolStoreKey(scope: String, name: String, iteration: Option[Int], input: Boolean)

object SymbolStoreEntry {
  def buildCallInputEntry(call: Call, inputName: String, wdlType: WdlType): SymbolStoreEntry = {
    val key = SymbolStoreKey(call.fullyQualifiedName, inputName, iteration = None, input = true)
    new SymbolStoreEntry(key, wdlType, None)
  }

  def buildCallOutputEntry(call: Call, taskOutput: TaskOutput): SymbolStoreEntry = {
    val key = SymbolStoreKey(call.fullyQualifiedName, taskOutput.name, iteration = None, input = false)
    new SymbolStoreEntry(key, taskOutput.wdlType, None)
  }
}

case class SymbolStoreEntry(key: SymbolStoreKey, wdlType: WdlType, wdlValue: Option[WdlValue]) {

  wdlValue match {
    case Some (value) => if(wdlType != value.wdlType) throw new UnsupportedOperationException("Fix this.")
    case None => // Nothing is always okay for a value.
  }

  def isInput: Boolean = key.input

  def isOutput: Boolean = !isInput

  /**
   * Create a copy of the receiver with `wdlValue` set to `newValue`.
   */
  def copyWithUpdatedValue(newValue: WdlValue): SymbolStoreEntry =
    new SymbolStoreEntry(key, wdlType, Option(newValue))

  def scope: String = key.scope
}

class SymbolStore(binding: WdlBinding, inputs: Map[FullyQualifiedName, WdlValue]) {
  private val store = mutable.Set[SymbolStoreEntry]()

  binding.workflow.calls.foreach { call =>
    call.task.inputs.foreach { case (inputName, wdlType) =>
      store += SymbolStoreEntry.buildCallInputEntry(call, inputName, wdlType)
    }
    call.task.outputs.foreach { output =>
      store += SymbolStoreEntry.buildCallOutputEntry(call, output)
    }
  }

  inputs.foreach { case (fullyQualifiedName, value) =>
    assignSymbolStoreEntry(fullyQualifiedName, value, iteration = None, input = true)
  }

  def locallyQualifiedInputs(call: Call): Map[FullyQualifiedName, WdlValue] = {
    store.collect { case entry if entry.scope == call.fullyQualifiedName && entry.isInput =>
      entry.key.name -> entry.wdlValue.get
    }.toMap
  }

  private def assignSymbolStoreEntry(fullyQualifiedName: FullyQualifiedName, value: WdlValue, iteration: Option[Int], input: Boolean): Unit = {
    val lastIndex = fullyQualifiedName.lastIndexOf(".")
    val (scope, name) = (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))

    val key = SymbolStoreKey(scope, name, iteration = None, input = true)

    val maybeOldEntry = store.find { _.key == key }

    maybeOldEntry match {
      case None =>
        throw new IllegalArgumentException(s"Failed to find expected symbol store entry with key '$key'")
      case Some(oldEntry) =>
        store.remove(oldEntry)
        store.add(oldEntry.copyWithUpdatedValue(value))
    }
  }
}
