package cwl

import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.Validation._
import cwl.CwlWomExpression.EnhancedParameterContextInputs
import wdl.types.{WdlMapType, WdlNothingType, WdlStringType, WdlType}
import wdl.values.{WdlFile, WdlMap, WdlString, WdlValue}
import wom.expression.{IoFunctionSet, WomExpression}

sealed trait CwlWomExpression extends WomExpression {

  def cwlExpressionType: WdlType

  override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = cwlExpressionType.validNel

  override def evaluateFiles(inputTypes: Map[String, WdlValue],
                             ioFunctionSet: IoFunctionSet,
                             coerceTo: WdlType): ErrorOr[Set[WdlFile]] = {
    ???
  }

  // TODO WOM oh geez
  override def inputs: Set[String] = ???
}

case class CommandOutputExpression(commandOutputParameter: CommandOutputParameter,
                                   override val cwlExpressionType: WdlType) extends CwlWomExpression {
  override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
    val parameterContext = ParameterContext.Empty.withInputs(inputValues, ioFunctionSet)

    commandOutputParameter.outputBinding match {
      case Some(outputBindingValue) =>
        val wdlValue = CommandOutputBindingEvaluator.commandOutputBindingToWdlValue(
          outputBindingValue,
          parameterContext,
          ioFunctionSet
        )
        cwlExpressionType.coerceRawValue(wdlValue).toErrorOr
      case None =>
        s"outputBinding not specified in $commandOutputParameter".invalidNel
    }
  }
}

case class StringExpression(expression: String) extends CwlWomExpression {
  override val cwlExpressionType = WdlStringType
  override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet) = {
    val parameterContext = ParameterContext.Empty.withInputs(inputValues, ioFunctionSet)
    // TODO: WOM: Instead of letting exceptions fly, catch and convert to ErrorOr
    ExpressionEvaluator.evalExpression(expression, parameterContext).valid
  }
}

object CwlWomExpression {
  def apply(commandOutputParameter: CommandOutputParameter, wdlType: WdlType): WomExpression = {
    CommandOutputExpression(commandOutputParameter, wdlType)
  }

  def apply(expression: String): WomExpression = {
    StringExpression(expression)
  }

  implicit class EnhancedParameterContextInputs(val parameterContext: ParameterContext) extends AnyVal {
    def withInputs(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ParameterContext = {
      val wdlValueType = inputValues.values.headOption.map(_.wdlType).getOrElse(WdlNothingType)
      parameterContext.copy(
        inputs = WdlMap(
          WdlMapType(WdlStringType, wdlValueType),
          // TODO: WOM: convert inputValues (including WdlFile?) to inputs using the ioFunctionSet
          inputValues map { case (name, wdlValue) => WdlString(name) -> wdlValue }
        )
      )
    }
  }
}
