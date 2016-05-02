package cromwell

import wdl4s.WdlExpression.ScopedLookupFunction
import wdl4s._
import wdl4s.expression.WdlFunctions
import wdl4s.values.{SymbolHash, WdlValue}

import scala.language.postfixOps
import scala.util.{Failure, Success, Try}
import scalaz._

package object core {
  // root can be a Path instead of a String in PBE. stdout / err too but it doesn't really bring values since they're just stringified to WdlFiles
  class WorkflowContext(val root: String)
  class CallContext(override val root: String, val stdout: String, val stderr: String) extends WorkflowContext(root)

  type ErrorOr[+A] = ValidationNel[String, A]
  type LocallyQualifiedName = String
  case class CallOutput(wdlValue: WdlValue, hash: Option[SymbolHash])
  type CallOutputs = Map[LocallyQualifiedName, CallOutput]
  type EvaluatedRuntimeAttributes = Map[String, WdlValue]

  implicit class EnhancedFQN(val fqn: FullyQualifiedName) extends AnyVal {
    def unqualified: LocallyQualifiedName = splitFqn._2

    def splitFqn: (String, String) = {
      val lastIndex = fqn.lastIndexOf(".")
      if (lastIndex == -1) ("", fqn)
      else (fqn.substring(0, lastIndex), fqn.substring(lastIndex + 1))
    }
  }
}