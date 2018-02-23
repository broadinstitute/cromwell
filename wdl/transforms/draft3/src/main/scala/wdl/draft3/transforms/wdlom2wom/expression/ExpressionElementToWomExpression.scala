package wdl.draft3.transforms.wdlom2wom.expression

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import wdl.draft3.transforms.wdlom2wom.expression
import wdl.draft3.transforms.wdlom2wom.linking.UnlinkedValueConsumer.ops._
import wdl.draft3.transforms.wdlom2wom.linking._
import wdl.draft3.transforms.wdlom2wom.expression.linking._
import wdl.model.draft3.elements.ExpressionElement
import wom.graph.GraphNodePort.OutputPort

object ExpressionElementToWomExpression {
  def convert(a: WomExpressionMakerInputs): ErrorOr[WdlomWomExpression] = {
    val consumedValues = a.expression.consumedValueNames

    def findConsumedPort(c: UnlinkedConsumedValueName): ErrorOr[LinkedConsumedValue] = c match {
      case c @ UnlinkedIdentifierName(sv) if a.availableInputs.contains(sv) =>
        LinkedConsumedValue(c, c, a.availableInputs(sv).womType).validNel
      case c @ UnlinkedCallOutputOrIdentifierAndMemberAccess(first, _) if a.availableInputs.contains(first) =>
        LinkedConsumedValue(c, UnlinkedIdentifierName(first), a.availableInputs(first).womType).validNel
      case c @ UnlinkedCallOutputOrIdentifierAndMemberAccess(first, second) if a.availableInputs.contains(s"$first.$second") =>
        LinkedConsumedValue(c, UnlinkedCallOutputName(first, second), a.availableInputs(s"$first.$second").womType).validNel
      case other => s"No lookup available for $other".invalidNel
    }

    val linkedConsumedValuesValidation: ErrorOr[Set[LinkedConsumedValue]] = consumedValues.toList traverse[ErrorOr, LinkedConsumedValue] { findConsumedPort } map (_.toSet)
    linkedConsumedValuesValidation.map(expression.WdlomWomExpression(a.expression, _))
  }
}

final case class WomExpressionMakerInputs(expression: ExpressionElement, availableInputs: Map[String, OutputPort])
