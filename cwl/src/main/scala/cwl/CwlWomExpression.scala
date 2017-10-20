package cwl

import cats.syntax.option._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cwl.WorkflowStepInput.InputSource
import shapeless.{Inl, Poly1}
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types._
import wom.values._

sealed trait CwlWomExpression extends WomExpression {

  def cwlExpressionType: WomType

  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = cwlExpressionType.validNel
}

case class CommandOutputExpression(outputBinding: CommandOutputBinding,
                                   override val cwlExpressionType: WomType,
                                   override val inputs: Set[String]) extends CwlWomExpression {

  // TODO WOM: outputBinding.toString is probably not be the best representation of the outputBinding
  override def sourceString = outputBinding.toString
  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
    val parameterContext = ParameterContext.Empty.withInputs(inputValues, ioFunctionSet)

    val wdlValue: WomValue = outputBinding.commandOutputBindingToWdlValue(parameterContext, ioFunctionSet)
    val extractFile: WomValue =
      wdlValue match {
        case WomArray(_, Seq(WomMap(WomMapType(WomStringType, WomStringType), map))) => map(WomString("location"))
        case other => other
      }
    cwlExpressionType.coerceRawValue(extractFile).toErrorOr
  }

  /*
  TODO:
   DB: It doesn't make sense to me that this function returns type WdlFile but accepts a type to which it coerces.
   Wouldn't coerceTo always == WdlFileType, and if not then what?
   */
  override def evaluateFiles(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] ={

    val pc = ParameterContext.Empty.withInputs(inputs, ioFunctionSet)

    val files = for {
      globValue <- outputBinding.glob.toList
      path <- GlobEvaluator.globPaths(globValue, pc, ioFunctionSet).toList
    } yield WomGlobFile(path): WomFile

    files.toSet.validNel[String]
  }
}

final case class WorkflowStepInputExpression(input: WorkflowStepInput, override val cwlExpressionType: WomType, graphInputs: Set[String]) extends CwlWomExpression {

  override def sourceString = input.toString

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) = {
    (input.valueFrom, input.source) match {
      case (None, Some(Inl(id: String))) =>
        inputValues.
          get(FullyQualifiedName(id).id).
          toValidNel(s"could not find id $id in typeMap\n${inputValues.mkString("\n")}\nwhen evaluating $input.  Graph Inputs were ${graphInputs.mkString("\n")}")
      case _ => s"Could not do evaluateValue(${input.valueFrom}, ${input.source}), most likely it has not been implemented yet".invalidNel
    }
  }

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) = ???

  object InputSourceToFileNames extends Poly1{

    implicit def string = at[String]{s => Set(FullyQualifiedName(s).id)}

    implicit def array = at[Array[String]]{_.map(FullyQualifiedName(_).id).toSet}
  }

  override def inputs = graphInputs ++ input.source.toSet.flatMap{(_:InputSource).fold(InputSourceToFileNames)}
}


