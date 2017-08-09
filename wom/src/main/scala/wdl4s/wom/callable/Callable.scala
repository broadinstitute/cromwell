package wdl4s.wom.callable

import lenthall.validation.ErrorOr.ErrorOr
import wdl4s.wdl.types.{WdlOptionalType, WdlType}
import wdl4s.wom.callable.Callable._
import wdl4s.wom.expression.WomExpression
import wdl4s.wom.graph.Graph


trait Callable {
  def name: String

  def graph: ErrorOr[Graph]
  def inputs: Set[_ <: InputDefinition]
  def declarations: List[(String, WomExpression)]
}

object Callable {
  sealed trait InputDefinition {
    def name: String
    def womType: WdlType
  }

  final case class OptionalInputDefinition(name: String, womType: WdlOptionalType) extends InputDefinition
  final case class OptionalInputDefinitionWithDefault(name: String, womType: WdlType, default: WomExpression) extends InputDefinition
  final case class RequiredInputDefinition(name: String, womType: WdlType) extends InputDefinition

  // Not really an input type, since it's not overridable by inputs but, shrug, I think it makes things easier to pretend they are inputs.
  final case class DeclaredInputDefinition(name: String, womType: WdlType, expression: WomExpression) extends InputDefinition

  case class OutputDefinition(name: String, womType: WdlType, expression: WomExpression)
}
