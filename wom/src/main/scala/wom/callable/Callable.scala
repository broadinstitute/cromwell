package wom.callable

import lenthall.validation.ErrorOr.ErrorOr
import wom.callable.Callable._
import wom.expression.WomExpression
import wom.graph.{Graph, LocalName}
import wom.types.{WdlOptionalType, WdlType}


trait Callable {
  def name: String

  def graph: ErrorOr[Graph]
  def inputs: List[_ <: InputDefinition]
}

object Callable {
  sealed trait InputDefinition {
    def localName: LocalName
    def womType: WdlType

    /**
      * Alias for localName.asString
      */
    def name = localName.value
  }

  object RequiredInputDefinition {
    def apply(name: String, womType: WdlType): RequiredInputDefinition = {
      RequiredInputDefinition(LocalName(name), womType)
    }
  }
  final case class RequiredInputDefinition(localName: LocalName, womType: WdlType) extends InputDefinition

  object InputDefinitionWithDefault {
    def apply(name: String, womType: WdlType, default: WomExpression): InputDefinitionWithDefault = {
      InputDefinitionWithDefault(LocalName(name), womType, default)
    }
  }
  final case class InputDefinitionWithDefault(localName: LocalName, womType: WdlType, default: WomExpression) extends InputDefinition

  object OptionalInputDefinition {
    def apply(name: String, womType: WdlOptionalType): OptionalInputDefinition = {
      OptionalInputDefinition(LocalName(name), womType)
    }
  }
  final case class OptionalInputDefinition(localName: LocalName, womType: WdlOptionalType) extends InputDefinition

  object OutputDefinition {
    def apply(name: String, womType: WdlType, expression: WomExpression): OutputDefinition = {
      OutputDefinition(LocalName(name), womType, expression)
    }
  }
  final case class OutputDefinition(localName: LocalName, womType: WdlType, expression: WomExpression) {
    /**
      * Alias for localName.asString
      */
    lazy val name = localName.value
  }
}
