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
}

case class CommandOutputExpression(outputBinding: CommandOutputBinding,
                                   override val cwlExpressionType: WdlType) extends CwlWomExpression {

  override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
    val parameterContext = ParameterContext.Empty.withInputs(inputValues, ioFunctionSet)

    outputBinding match {
      case outputBindingValue =>
        val wdlValue = CommandOutputBindingEvaluator.commandOutputBindingToWdlValue(
          outputBindingValue,
          parameterContext,
          ioFunctionSet
        )
        cwlExpressionType.coerceRawValue(wdlValue).toErrorOr
    }
  }

  override def inputs: Set[String] = ???

  override def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] = ???
}

object CwlWomExpression {

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
