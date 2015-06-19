package cromwell.engine.store

import cromwell.binding._
import cromwell.binding.types.WdlType
import cromwell.binding.values.{WdlObject, WdlValue}
import cromwell.engine.store.SymbolStore._

import scala.util.{Failure, Success, Try}

object SymbolStore {
  case class SymbolStoreKey(scope: String, name: String, iteration: Option[Int], input: Boolean)

  object CallInputWdlFunctions extends WdlFunctions {
    def getFunction(name: String): WdlFunction = {
      throw new WdlExpressionException("TODO: Some functions may be allowed in this context")
    }
  }

  def apply(namespace: WdlNamespace, inputs: HostInputs): SymbolStore = {
    new SymbolStore(namespace, inputs, Set.empty[SymbolStoreEntry])
  }
}

case class SymbolStore(namespace: WdlNamespace, inputs: HostInputs, initialStore: Set[SymbolStoreEntry]) {
  private var store: Set[SymbolStoreEntry] = createStore()

  def locallyQualifiedInputs(call: Call): Map[String, WdlValue] = {
    def lookup(identifierString: String): WdlValue = {
      val workflow = call.parent.map {_.asInstanceOf[Workflow]} getOrElse {
        throw new WdlExpressionException("Expecting 'call' to have a 'workflow' parent")
      }
      val namespaces = call.namespace.namespaces filter {_.namespace.contains(identifierString)}
      namespaces.headOption.getOrElse {
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
        case Some(e: WdlExpression) => e.evaluate(lookup, CallInputWdlFunctions).get
        case Some(v) => v
        case _ => throw new WdlExpressionException("Unknown error")
      }
      entry.key.name -> value
    }.toMap
  }

  def outputs = store filter {_.isOutput}

  def addOutputValue(scope: String, name: String, maybeValue: Option[WdlValue], wdlType: WdlType): Try[Unit] = {
    val key = SymbolStoreKey(scope, name, iteration = None, input = false)
    store find {_.key == key} match {
      case Some(preexisting) =>
        Failure(new IllegalArgumentException(s"Found unexpected preexisting symbol store entry with key '$key'"))
      case None => Success(store += SymbolStoreEntry(key, wdlType, maybeValue)) // Success of a side effect? FTW!
    }
  }

  private def callEntries(call: Call) = store.filter {_.scope == call.fullyQualifiedName}
  private def callOutputEntries(call: Call) = callEntries(call).filter {entry => entry.isOutput}
  private def callInputEntries(call: Call) = callEntries(call).filter {entry => entry.isInput}

  private def createStore(): Set[SymbolStoreEntry] = if (initialStore.nonEmpty) initialStore else createInitialStore

  private def createInitialStore: Set[SymbolStoreEntry] = {
    // FIXME: Is this terminology correct?
    val inputSymbols = inputs.map {case (name, value) => SymbolStoreEntry(name, value, input = true)}
    val callSymbols = for {
      workflow <- namespace.workflows
      call <- workflow.calls
      (k, v) <- call.inputMappings
    } yield SymbolStoreEntry(s"${call.fullyQualifiedName}.$k", v, input = true)

    inputSymbols.toSet ++ callSymbols.toSet
  }

  override def toString: String = {
    store.map{e => s"${e.key.scope}\t${e.key.name}\t${e.key.input}\t${e.wdlType}\t${e.wdlValue}"}.mkString("\n")
  }
}


