package cwl

import cats.syntax.validated._
import lenthall.validation.ErrorOr.ErrorOr
import lenthall.validation.Validation._
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types._
import wom.values._

sealed trait CwlWomExpression extends WomExpression {

  def cwlExpressionType: WdlType

  override def evaluateType(inputTypes: Map[String, WdlType]): ErrorOr[WdlType] = cwlExpressionType.validNel
}

case class CommandOutputExpression(outputBinding: CommandOutputBinding,
                                   override val cwlExpressionType: WdlType) extends CwlWomExpression {

  // TODO WOM: outputBinding.toString is probably not be the best representation of the outputBinding
  override def sourceString = outputBinding.toString
  override def evaluateValue(inputValues: Map[String, WdlValue], ioFunctionSet: IoFunctionSet): ErrorOr[WdlValue] = {
    val parameterContext = ParameterContext.Empty.withInputs(inputValues, ioFunctionSet)

    val wdlValue = outputBinding.commandOutputBindingToWdlValue(parameterContext, ioFunctionSet)
    cwlExpressionType.coerceRawValue(wdlValue).toErrorOr
  }

  override def inputs: Set[String] = ???

  /*
  TODO:
   DB: It doesn't make sense to me that this function returns type WdlFile but accepts a type to which it coerces.
   Wouldn't coerceTo always == WdlFileType, and if not then what?
   */
  override def evaluateFiles(inputTypes: Map[String, WdlValue], ioFunctionSet: IoFunctionSet, coerceTo: WdlType): ErrorOr[Set[WdlFile]] ={

    val pc = ParameterContext.Empty.withInputs(inputTypes, ioFunctionSet)
    val wdlValue = outputBinding.commandOutputBindingToWdlValue(pc, ioFunctionSet)

    wdlValue match {

      case WdlArray(WdlMaybeEmptyArrayType(WdlMapType(WdlStringType, WdlStringType)), seq: Seq[WdlValue]) =>
        seq.map {
          case WdlMap(WdlMapType(WdlStringType, WdlStringType), map) => WdlGlobFile(map(WdlString("location")).valueString): WdlFile
        }.toSet.validNel

      case other =>s":( we saw $other and couldn't convert to a globfile type: ${other.wdlType} coerceTo: $coerceTo".invalidNel[Set[WdlFile]]
    }
  }
}
