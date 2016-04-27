package cromwell

import wdl4s._
import wdl4s.expression.WdlFunctions
import wdl4s.values.{SymbolHash, WdlValue}

import scala.util.Try
import scalaz._

package object core {
  // root can be a Path instead of a String in PBE. stdout / err too but it doesn't really bring values since they're just stringified to WdlFiles
  class WorkflowContext(val root: String)
  class CallContext(override val root: String, val stdout: String, val stderr: String) extends WorkflowContext(root)

  // Can evaluate a wdl value
  class Evaluator(evaluator: WdlValue => Try[WdlValue]) {
    def evaluate(wdlValue: WdlValue): Try[WdlValue] = evaluator(wdlValue)
  }

  // Can build an evaluator from engine functions and valueMapper
  class EvaluatorBuilder(builder: (WdlFunctions[WdlValue], WdlValue => WdlValue) => WdlValue => Try[WdlValue]) {
    def build(engineFunctions: WdlFunctions[WdlValue], preValueMapper: WdlValue => WdlValue = identity): Evaluator = {
      new Evaluator(builder(engineFunctions, preValueMapper))
    }
  }

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