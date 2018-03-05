package cwl

import cats.instances.option._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import shapeless.{:+:, CNil, Poly1}
import wom.values.{WomString, WomValue}

trait OutputParameter {
  def id: String
  def label: Option[String]
  def secondaryFiles: Option[SecondaryFiles]
  def format: Option[OutputParameterFormat]
  def streamable: Option[Boolean]
  def doc: Option[String :+: Array[String] :+: CNil]
  def outputBinding: Option[CommandOutputBinding]
  def `type`: Option[MyriadOutputType]
}

object OutputParameter {
  object IdAndType {
    def unapply(arg: OutputParameter): Option[(String, MyriadOutputType)] = arg.`type`.map((arg.id, _))
  }

  def format(formatOption: Option[StringOrExpression],
             parameterContext: ParameterContext,
             expressionLib: ExpressionLib): ErrorOr[Option[String]] = {
    formatOption.traverse[ErrorOr, String] {
      format(_, parameterContext, expressionLib)
    }
  }

  def format(format: StringOrExpression, parameterContext: ParameterContext, expressionLib: ExpressionLib): ErrorOr[String] = {
    format.fold(OutputParameter.FormatPoly).apply(parameterContext, expressionLib)
  }

  type FormatFunction = (ParameterContext, ExpressionLib) => ErrorOr[String]

  object FormatPoly extends Poly1 {
    implicit def caseStringOrExpression: Case.Aux[StringOrExpression, FormatFunction] = {
      at {
        _.fold(this)
      }
    }

    implicit def caseExpression: Case.Aux[Expression, FormatFunction] = {
      at {
        expression =>
          (parameterContext, expressionLib) =>
            val result: ErrorOr[WomValue] = ExpressionEvaluator.eval(expression, parameterContext, expressionLib)
            result flatMap {
              case womString: WomString => womString.value.valid
              case other => s"Not a valid file format: $other".invalidNel
            }
      }
    }

    implicit def caseString: Case.Aux[String, FormatFunction] = at { string => (_,_) => string.valid }
  }
}
