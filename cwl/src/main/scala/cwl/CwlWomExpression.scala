package cwl

import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types._
import wom.values._

sealed trait CwlWomExpression extends WomExpression {

  def cwlExpressionType: WomType

  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = cwlExpressionType.validNel
}

case class CommandOutputExpression(outputBinding: CommandOutputBinding,
                                   override val cwlExpressionType: WomType) extends CwlWomExpression {

  // TODO WOM: outputBinding.toString is probably not be the best representation of the outputBinding
  override def sourceString = outputBinding.toString
  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
    val parameterContext = ParameterContext.Empty.withInputs(inputValues, ioFunctionSet)

    val womValue = outputBinding.commandOutputBindingToWdlValue(parameterContext, ioFunctionSet)
    cwlExpressionType.coerceRawValue(womValue).toErrorOr
  }

  override def inputs: Set[String] = ???

  /*
  TODO:
   DB: It doesn't make sense to me that this function returns type WdlFile but accepts a type to which it coerces.
   Wouldn't coerceTo always == WdlFileType, and if not then what?
   */
  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] ={

    val pc = ParameterContext.Empty.withInputs(inputTypes, ioFunctionSet)
    val womValue = outputBinding.commandOutputBindingToWdlValue(pc, ioFunctionSet)

    womValue match {

      case WomArray(WomMaybeEmptyArrayType(WomMapType(WomStringType, WomStringType)), seq: Seq[WomValue]) =>
        seq.map {
          case WomMap(WomMapType(WomStringType, WomStringType), map) => WomGlobFile(map(WomString("location")).valueString): WomFile
        }.toSet.validNel

      case other =>s":( we saw $other and couldn't convert to a globfile type: ${other.womType} coerceTo: $coerceTo".invalidNel[Set[WomFile]]
    }
  }
}
