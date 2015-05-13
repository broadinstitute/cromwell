package cromwell.engine

import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.WdlValue

import scala.collection.mutable
import scala.util.{Success, Failure, Try}


case class SymbolStoreKey(scope: String, name: String, iteration: Option[Int], input: Boolean)

object SymbolStoreEntry {
  def buildCallInputEntry(call: Call, inputName: String, wdlType: WdlType): SymbolStoreEntry = {
    val key = SymbolStoreKey(call.fullyQualifiedName, inputName, iteration = None, input = true)
    SymbolStoreEntry(key, wdlType, None)
  }

  def buildCallOutputEntry(call: Call, taskOutput: TaskOutput): SymbolStoreEntry = {
    val key = SymbolStoreKey(call.fullyQualifiedName, taskOutput.name, iteration = None, input = false)
    SymbolStoreEntry(key, taskOutput.wdlType, None)
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
  }

  inputs.foreach { case (fullyQualifiedName, value) =>
    assignSymbolStoreEntry(fullyQualifiedName, value, iteration = None, input = true)
  }

  def locallyQualifiedInputs(call: Call): Map[String, WdlValue] = {
    store.collect { case entry if entry.scope == call.fullyQualifiedName && entry.isInput =>
      entry.key.name -> entry.wdlValue.get
    }.toMap
  }

  private def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
    val lastIndex = fullyQualifiedName.lastIndexOf(".")
    (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
  }

  private def assignSymbolStoreEntry(fullyQualifiedName: FullyQualifiedName, value: WdlValue, iteration: Option[Int], input: Boolean): Unit = {
    val (scope, name) = splitFqn(fullyQualifiedName)

    val key = SymbolStoreKey(scope, name, iteration = None, input = true)

    store.find { _.key == key} match {
      case None =>
        throw new IllegalArgumentException(s"Failed to find expected symbol store entry with key '$key'")
      case Some(oldEntry) =>
        store.remove(oldEntry)
        store.add(oldEntry.copyWithUpdatedValue(value))
    }
  }
  
  def addOutputValue(scope: String, name: String, maybeValue: Option[WdlValue], wdlType: WdlType): Try[Unit] = {
    val key = SymbolStoreKey(scope, name, iteration = None, input = false)
    
    val maybeOldEntry = store.find { _.key == key }
    maybeOldEntry match {
      case Some(preexisting) =>
        Failure(new IllegalArgumentException(s"Found unexpected preexisting symbol store entry with key '$key'"))
      case None =>
        Success(store.add(SymbolStoreEntry(key, wdlType, maybeValue)))
    }
  }

  def getOutputByFullyQualifiedName(fqn: FullyQualifiedName): Option[SymbolStoreEntry] = {
    val (scope, name) = splitFqn(fqn)
    val key = SymbolStoreKey(scope, name, iteration = None, input = false)
    store.find(_.key == key)
  }
}
