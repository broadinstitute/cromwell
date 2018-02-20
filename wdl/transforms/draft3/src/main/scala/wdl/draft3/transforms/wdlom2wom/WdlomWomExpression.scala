package wdl.draft3.transforms.wdlom2wom

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.wdlom2wom.ValueConsumer.ops._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{IdentifierLookup, PrimitiveLiteralExpressionElement}
import wom.expression.{IoFunctionSet, WomExpression}
import wom.graph.GraphNodePort.OutputPort
import wom.types.WomType
import wom.values.{WomFile, WomValue}

final case class WdlomWomExpression(expressionElement: ExpressionElement, inputMap: Map[ConsumedValue, FoundConsumedValue]) extends WomExpression {
  override def sourceString: String = expressionElement.toString

  override def inputs: Set[String] = inputMap.values.map(_.inputName).toSet

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = expressionElement match {
    case PrimitiveLiteralExpressionElement(v) => v.validNel
    case IdentifierLookup(identifier) => inputValues.get(identifier).map(_.validNel).getOrElse(s"No value supplied for variable lookup: $identifier".invalidNel)
    case _ => ??? // TODO other expression elements
  }

  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = expressionElement match {
    case PrimitiveLiteralExpressionElement(v) => v.womType.validNel
    case id: IdentifierLookup => inputMap((id: ExpressionElement).consumedValues.head).outputPort.womType.validNel
    case _ => ??? // TODO other expression elements
  }

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] = {
    Set.empty[WomFile].validNel
  }
}

final case class FoundConsumedValue(outputPort: OutputPort, inputName: String)
