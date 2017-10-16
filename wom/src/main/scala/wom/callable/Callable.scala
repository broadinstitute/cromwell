package wom.callable

import lenthall.validation.ErrorOr.ErrorOr
import wom.callable.Callable._
import wom.expression.WomExpression
import wom.graph.{Graph, LocalName}
import wom.types.{WomOptionalType, WomType}


trait Callable {
  def name: String

  def graph: ErrorOr[Graph]
  def inputs: List[_ <: InputDefinition]
}

object Callable {
  sealed trait InputDefinition {
    def localName: LocalName
    def womType: WomType

    /**
      * Alias for localName.asString
      */
    def name = localName.value
  }

  object RequiredInputDefinition {
    def apply(name: String, womType: WomType): RequiredInputDefinition = {
      RequiredInputDefinition(LocalName(name), womType)
    }
  }
  final case class RequiredInputDefinition(localName: LocalName, womType: WomType) extends InputDefinition

  object InputDefinitionWithDefault {
    def apply(name: String, womType: WomType, default: WomExpression): InputDefinitionWithDefault = {
      InputDefinitionWithDefault(LocalName(name), womType, default)
    }
  }
  final case class InputDefinitionWithDefault(localName: LocalName, womType: WomType, default: WomExpression) extends InputDefinition

  object OptionalInputDefinition {
    def apply(name: String, womType: WomOptionalType): OptionalInputDefinition = {
      OptionalInputDefinition(LocalName(name), womType)
    }
  }
  final case class OptionalInputDefinition(localName: LocalName, womType: WomOptionalType) extends InputDefinition

  object OutputDefinition {
    def apply(name: String, womType: WomType, expression: WomExpression): OutputDefinition = {
      OutputDefinition(LocalName(name), womType, expression)
    }
  }
  final case class OutputDefinition(localName: LocalName, womType: WomType, expression: WomExpression) {
    /**
      * Alias for localName.asString
      */
    lazy val name = localName.value
  }
}
