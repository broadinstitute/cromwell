package cromwell.engine

import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.WdlValue

import scala.collection.mutable
import scala.util.{Failure, Success, Try}


case class SymbolStoreKey(scope: String, name: String, iteration: Option[Int], input: Boolean)

case class SymbolStoreEntry(key: SymbolStoreKey, wdlType: WdlType, wdlValue: Option[WdlValue]) {

  def isInput: Boolean = key.input

  def isOutput: Boolean = !isInput

  def scope: String = key.scope
}

class SymbolStore(binding: WdlBinding, inputs: Map[FullyQualifiedName, WdlValue]) {
  private val store = mutable.Set[SymbolStoreEntry]()

  inputs.foreach { case (fullyQualifiedName, value) =>
    assignSymbolStoreEntry(fullyQualifiedName, value, input = true)
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

  private def assignSymbolStoreEntry(fullyQualifiedName: FullyQualifiedName, wdlValue: WdlValue, input: Boolean): Unit = {
    val (scope, name) = splitFqn(fullyQualifiedName)

    val key = SymbolStoreKey(scope, name, iteration = None, input = true)
    store.add(SymbolStoreEntry(key, wdlValue.wdlType, Some(wdlValue)))
  }

  def addOutputValue(scope: String, name: String, maybeValue: Option[WdlValue], wdlType: WdlType): Try[Unit] = {
    val key = SymbolStoreKey(scope, name, iteration = None, input = false)

    val maybeOldEntry = store.find {
      _.key == key
    }
    maybeOldEntry match {
      case Some(preexisting) =>
        Failure(new IllegalArgumentException(s"Found unexpected preexisting symbol store entry with key '$key'"))
      case None =>
        Success(store.add(SymbolStoreEntry(key, wdlType, maybeValue)))
    }
  }

  def addInputValue(scope: String, name: String, entry: SymbolStoreEntry): Unit = {
    val key = SymbolStoreKey(scope, name, iteration = None, input = true)
    store.add(SymbolStoreEntry(key, entry.wdlType, entry.wdlValue))
  }

  def copyOutputToInput(outputFqn: String, inputFqn: String): Try[Unit] = {
    getOutputByFullyQualifiedName(outputFqn) match {
      case None =>
        Failure(new IllegalArgumentException(s"Did not find expected output '$outputFqn'"))
      case Some(entry) =>
        val (inputScope, inputName) = splitFqn(inputFqn)
        Success(addInputValue(inputScope, inputName, entry))
    }
  }

  def getOutputByFullyQualifiedName(fqn: FullyQualifiedName): Option[SymbolStoreEntry] = {
    val (scope, name) = splitFqn(fqn)
    val key = SymbolStoreKey(scope, name, iteration = None, input = false)
    store.find(_.key == key)
  }

  def getOutputs: Set[SymbolStoreEntry] = store.toSet filter {_.isOutput}
}
