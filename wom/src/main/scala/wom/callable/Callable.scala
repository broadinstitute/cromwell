package wom.callable

import common.validation.ErrorOr.ErrorOr
import wom.callable.Callable.InputDefinition.InputValueMapper
import wom.callable.Callable._
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.{CommandCallNode, Graph, LocalName}
import wom.types.{WomOptionalType, WomType}
import wom.values.WomValue

trait Callable {
  def name: String

  def inputs: List[_ <: InputDefinition]
  def outputs: List[_ <: OutputDefinition]
}

trait ExecutableCallable extends Callable {
  def graph: Graph
  lazy val taskCallNodes: Set[CommandCallNode] = graph.allNodes collect {
    case taskNode: CommandCallNode => taskNode
  }
}

object Callable {
  object InputDefinition {
    type InputValueMapper = IoFunctionSet => WomValue => WomValue
    val IdentityValueMapper: InputValueMapper = { _ => value => value}
  }
  sealed trait InputDefinition {
    def localName: LocalName
    def womType: WomType

    /**
      * Alias for localName.asString
      */
    def name = localName.value
    def optional = this match {
      case _: RequiredInputDefinition => false
      case _ => true
    }
    def valueMapper: InputValueMapper
  }

  object RequiredInputDefinition {
    def apply(name: String, womType: WomType, valueMapper: InputValueMapper): RequiredInputDefinition = {
      RequiredInputDefinition(LocalName(name), womType, valueMapper)
    }
    def apply(name: String, womType: WomType): RequiredInputDefinition = {
      RequiredInputDefinition(LocalName(name), womType, InputDefinition.IdentityValueMapper)
    }
  }
  final case class RequiredInputDefinition(localName: LocalName, womType: WomType, valueMapper: InputValueMapper = InputDefinition.IdentityValueMapper) extends InputDefinition

  object InputDefinitionWithDefault {
    def apply(name: String, womType: WomType, default: WomExpression): InputDefinitionWithDefault = {
      InputDefinitionWithDefault(LocalName(name), womType, default, InputDefinition.IdentityValueMapper)
    }
    def apply(name: String, womType: WomType, default: WomExpression, valueMapper: InputValueMapper): InputDefinitionWithDefault = {
      InputDefinitionWithDefault(LocalName(name), womType, default, valueMapper)
    }
  }

  /**
    * An input definition that has a default value supplied. Typical WDL example would be a declaration like: "Int x = 5"
    */
  final case class InputDefinitionWithDefault(localName: LocalName, womType: WomType, default: WomExpression, valueMapper: InputValueMapper = InputDefinition.IdentityValueMapper) extends InputDefinition

  object OptionalInputDefinition {
    def apply(name: String, womType: WomOptionalType): OptionalInputDefinition = OptionalInputDefinition(LocalName(name), womType, InputDefinition.IdentityValueMapper)
    def apply(name: String, womType: WomOptionalType, valueMapper: InputValueMapper): OptionalInputDefinition = OptionalInputDefinition(LocalName(name), womType, valueMapper)
  }
  final case class OptionalInputDefinition(localName: LocalName, womType: WomOptionalType, valueMapper: InputValueMapper = InputDefinition.IdentityValueMapper) extends InputDefinition

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
