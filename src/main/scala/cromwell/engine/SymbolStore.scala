package cromwell.engine

import scala.collection.{mutable, breakOut}
import scala.util.{Failure, Success, Try}
import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.{WdlObject, WdlValue}
import SymbolStore._

object SymbolStore {
  case class SymbolStoreKey(scope: String, name: String, iteration: Option[Int], input: Boolean)

  case class SymbolStoreEntry(key: SymbolStoreKey, wdlType: WdlType, wdlValue: Option[WdlValue]) {
    def isInput: Boolean = key.input
    def isOutput: Boolean = !isInput
    def scope: String = key.scope
  }

  object CallInputWdlFunctions extends WdlFunctions {
    def getFunction(name: String): WdlFunction = {
      throw new WdlExpressionException("TODO: Some functions may be allowed in this context")
    }
  }

  def apply(namespace: WdlNamespace, inputs: HostInputs): SymbolStore = new SymbolStore(namespace, inputs, Set.empty[SymbolStoreEntry])
}

case class SymbolStore(namespace: WdlNamespace, inputs: HostInputs, initialStore: Set[SymbolStoreEntry]) {
  // Convert the incoming Seq to a mutable Set
  private val store: mutable.Set[SymbolStoreEntry] = initialStore.map(identity)(breakOut)

  inputs.foreach { case (fullyQualifiedName, value) =>
    assignSymbolStoreEntry(fullyQualifiedName, value, input = true)
  }

  for {
    workflow <- namespace.workflows
    call <- workflow.calls
    (k, v) <- call.inputMappings
  } yield assignSymbolStoreEntry(s"${call.fullyQualifiedName}.$k", v, input = true)

  private def callEntries(call: Call): Set[SymbolStoreEntry] =
    store.filter {_.scope == call.fullyQualifiedName}.toSet

  private def callOutputEntries(call: Call): Set[SymbolStoreEntry] =
    callEntries(call).filter {entry => entry.isOutput}

  private def callInputEntries(call: Call): Set[SymbolStoreEntry] =
    callEntries(call).filter {entry => entry.isInput}

  def locallyQualifiedInputs(call: Call): Map[String, WdlValue] = {
    def lookup(identifierString: String): WdlValue = {
      val workflow = call.parent.map {_.asInstanceOf[Workflow]} getOrElse {
        throw new WdlExpressionException("Expecting 'call' to have a 'workflow' parent")
      }
      val namespaces = call.namespace.namespaces.filter{_.namespace.contains(identifierString)}
      namespaces.headOption.getOrElse{
        val matchedCall = workflow.calls.find {_.name == identifierString}.getOrElse {
          throw new WdlExpressionException(s"Expecting to find a call with name '$identifierString'")
        }
        val callOutputs = callOutputEntries(matchedCall) map { entry =>
          val value = entry.wdlValue match {
            case Some(v) => v
            case _ => throw new WdlExpressionException(s"Could not evaluate call '${matchedCall.name}', because '${entry.key.name}' is undefined")
          }
          entry.key.name -> value
        }
        WdlObject(callOutputs.toMap)
      }
    }

    callInputEntries(call).map {entry =>
      val value = entry.wdlValue match {
        case Some(e: WdlExpression) => e.evaluate(lookup, SymbolStore.CallInputWdlFunctions).get
        case Some(v) => v
        case _ => throw new WdlExpressionException("Unknown error")
      }
      entry.key.name -> value
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

  def getOutputs: Set[SymbolStoreEntry] = store.toSet filter {_.isOutput}

  def readOnly = store.toSet

  override def toString: String = {
    store.map{e => s"${e.key.scope}\t${e.key.name}\t${e.key.input}\t${e.wdlType}\t${e.wdlValue}"}.mkString("\n")
  }
}


