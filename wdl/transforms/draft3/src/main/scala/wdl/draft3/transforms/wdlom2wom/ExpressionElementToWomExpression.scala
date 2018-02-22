package wdl.draft3.transforms.wdlom2wom

import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._

import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wom.graph.GraphNodePort.OutputPort
import wdl.draft3.transforms.wdlom2wom.ValueConsumer.ops._

object ExpressionElementToWomExpression {
  def convert(a: WomExpressionMakerInputs): ErrorOr[WdlomWomExpression] = {
    val consumedValues = a.expression.consumedValues

    def findConsumedPort(c: ConsumedValue): ErrorOr[(ConsumedValue, FoundConsumedValue)] = c match {
      case c @ ConsumedSingleValue(sv) if a.availableInputs.contains(sv) => (c -> FoundConsumedValue(sv, a.availableInputs(sv).womType)).validNel
      case c @ ConsumedLookupValue(first, _) if a.availableInputs.contains(first) => (c -> FoundConsumedValue(first, a.availableInputs(first).womType)).validNel
      case c @ ConsumedLookupValue(first, second) if a.availableInputs.contains(s"$first.$second") => (c -> FoundConsumedValue(s"$first.$second", a.availableInputs(s"$first.$second").womType)).validNel
      case other => s"No lookup available for $other".invalidNel
    }

    val consumptionMapValidation: ErrorOr[Map[ConsumedValue, FoundConsumedValue]] = (consumedValues.toList traverse[ErrorOr, (ConsumedValue, FoundConsumedValue)] { findConsumedPort }).map(_.toMap)
    consumptionMapValidation.map(WdlomWomExpression(a.expression, _))
  }
}

final case class WomExpressionMakerInputs(expression: ExpressionElement, availableInputs: Map[String, OutputPort])
