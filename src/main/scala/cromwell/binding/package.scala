package cromwell

import cromwell.binding.values.{WdlFile, WdlValue}

import scala.language.implicitConversions

/**
 * ==WDL Bindings for Scala==
 *
 * This package contains an API to convert WDL ASTs into native Scala objects.
 *
 * Clients of this package should be able to interface with a WDL file without ever coming in direct contact with an Abstract Syntax Tree
 *
 * Internally, this package is built on top of [[cromwell.parser]].
 */

package object binding {
  type WdlSource = String
  type WdlJson = String
  type WorkflowOptionsJson = String
  type WorkflowRawInputs = Map[FullyQualifiedName, Any]
  type WorkflowCoercedInputs = Map[FullyQualifiedName, WdlValue]
  type WorkflowOutputs = Map[FullyQualifiedName, CallOutput]
  type FullyQualifiedName = String
  type LocallyQualifiedName = String
  type CallInputs = Map[String, WdlValue]
  case class CallOutput(wdlValue: WdlValue, hash: SymbolHash)
  type CallOutputs = Map[LocallyQualifiedName, CallOutput]
  type HostInputs = Map[String, WdlValue]
  type ImportResolver = String => WdlSource
  case class SymbolHash(value: String) extends Ordered[SymbolHash] {
    def compare(that: SymbolHash) = this.value compare that.value
  }
  type FileHasher = WdlFile => SymbolHash

  implicit class EnhancedFullyQualifiedName(val fqn: FullyQualifiedName) extends AnyVal {
    def isScatter = fqn.contains(Scatter.FQNIdentifier)
    def scopeAndVariableName: (String, String) = {
      val array = fqn.split("\\.(?=[^\\.]+$)")
      (array(0), array(1))
    }
  }

  implicit class EnhancedCallOutputMap[A](val m: Map[A, CallOutput]) extends AnyVal {
    def mapToValues: Map[A, WdlValue] = m map {
      case (k, CallOutput(wdlValue, hash)) => (k, wdlValue)
    }
  }
}
