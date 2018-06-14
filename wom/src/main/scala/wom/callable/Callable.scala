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
    import cats.syntax.validated._
    type InputValueMapper = IoFunctionSet => WomValue => ErrorOr[WomValue]
    val IdentityValueMapper: InputValueMapper = { _ => value => value.validNel}
  }
  sealed trait InputDefinition {
    def localName: LocalName
    def womType: WomType
    def parameterMeta: Option[MetaValueElement]

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
      RequiredInputDefinition(LocalName(name), womType, valueMapper, None)
    }
    def apply(name: String, womType: WomType): RequiredInputDefinition = {
      RequiredInputDefinition(LocalName(name), womType, InputDefinition.IdentityValueMapper, None)
    }
    def apply(name: String, womType: WomType, parameterMeta: Option[MetaValueElement]): RequiredInputDefinition = {
      RequiredInputDefinition(LocalName(name), womType, InputDefinition.IdentityValueMapper, parameterMeta)
    }
  }
  final case class RequiredInputDefinition(localName: LocalName,
                                           womType: WomType,
                                           valueMapper: InputValueMapper = InputDefinition.IdentityValueMapper,
                                           parameterMeta: Option[MetaValueElement] = None) extends InputDefinition

  object InputDefinitionWithDefault {
    def apply(name: String, womType: WomType, default: WomExpression): InputDefinitionWithDefault = {
      InputDefinitionWithDefault(LocalName(name), womType, default, InputDefinition.IdentityValueMapper, None)
    }
    def apply(name: String, womType: WomType, default: WomExpression, parameterMeta: Option[MetaValueElement]): InputDefinitionWithDefault = {
      InputDefinitionWithDefault(LocalName(name), womType, default, InputDefinition.IdentityValueMapper, parameterMeta)
    }
    def apply(name: String, womType: WomType, default: WomExpression, valueMapper: InputValueMapper): InputDefinitionWithDefault = {
      InputDefinitionWithDefault(LocalName(name), womType, default, valueMapper, None)
    }
  }

  /**
    * An input definition that has a default value supplied. Typical WDL example would be a declaration like: "Int x = 5"
    */
  final case class InputDefinitionWithDefault(localName: LocalName,
                                              womType: WomType,
                                              default: WomExpression,
                                              valueMapper: InputValueMapper = InputDefinition.IdentityValueMapper,
                                              parameterMeta: Option[MetaValueElement] = None) extends InputDefinition

  object FixedInputDefinition {
    def apply(name: String, womType: WomType, default: WomExpression): FixedInputDefinition = {
      FixedInputDefinition(LocalName(name), womType, default, InputDefinition.IdentityValueMapper, None)
    }
    def apply(name: String, womType: WomType, default: WomExpression, parameterMeta: Option[MetaValueElement]): FixedInputDefinition = {
      FixedInputDefinition(LocalName(name), womType, default, InputDefinition.IdentityValueMapper, parameterMeta)
    }
  }

  /**
    * An input whose value should always be calculated from the default, and is not allowed to be overridden.
    */
  final case class FixedInputDefinition(localName: LocalName,
                                        womType: WomType,
                                        default: WomExpression,
                                        valueMapper: InputValueMapper = InputDefinition.IdentityValueMapper,
                                        parameterMeta: Option[MetaValueElement] = None) extends InputDefinition

  object OptionalInputDefinition {
    def apply(name: String, womType: WomOptionalType): OptionalInputDefinition = OptionalInputDefinition(LocalName(name), womType, InputDefinition.IdentityValueMapper, None)
    def apply(name: String, womType: WomOptionalType, parameterMeta: Option[MetaValueElement]): OptionalInputDefinition = OptionalInputDefinition(LocalName(name), womType, InputDefinition.IdentityValueMapper, parameterMeta)
    def apply(name: String, womType: WomOptionalType, valueMapper: InputValueMapper): OptionalInputDefinition = OptionalInputDefinition(LocalName(name), womType, valueMapper, None)
  }
  final case class OptionalInputDefinition(localName: LocalName,
                                           womType: WomOptionalType,
                                           valueMapper: InputValueMapper = InputDefinition.IdentityValueMapper,
                                           parameterMeta: Option[MetaValueElement] = None) extends InputDefinition

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
