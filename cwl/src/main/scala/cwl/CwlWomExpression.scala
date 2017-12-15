package cwl

import cats.syntax.option._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wom.types._
import wom.values._
import wom.expression.{IoFunctionSet, WomExpression}
import cats.syntax.validated._
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName

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

    //To facilitate ECMAScript evaluation, filenames are stored in a map under the key "location"
    val womValue = outputBinding.
      commandOutputBindingToWomValue(parameterContext, ioFunctionSet) match {
        case WomArray(_, Seq(WomMap(WomMapType(WomStringType, WomStringType), map))) => map(WomString("location"))
        case other => other
      }

    //If the value is a string but the output is expecting a file, we consider that string a POSIX "glob" and apply
    //it accordingly to retrieve the file list to which it expands.
    val globbedIfFile =
      (womValue, cwlExpressionType) match {

        //In the case of a single file being expected, we must enforce that the glob only represents a single file
        case (WomString(glob), WomSingleFileType) =>
          ioFunctionSet.glob(glob) match {
            case head :: Nil => WomString(head)
            case list => throw new RuntimeException(s"expecting a single File glob but instead got $list")
          }

        case _  => womValue
      }

    //CWL tells us the type this output is expected to be.  Attempt to coerce the actual output into this type.
    cwlExpressionType.coerceRawValue(globbedIfFile).toErrorOr
  }

  /*
  TODO:
   DB: It doesn't make sense to me that this function returns type WomFile but accepts a type to which it coerces.
   Wouldn't coerceTo always == WomFileType, and if not then what?
   */
  override def evaluateFiles(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] ={

    val pc = ParameterContext().withInputs(inputs, ioFunctionSet)

    val files = for {
      globValue <- outputBinding.glob.toList
      path <- GlobEvaluator.globPaths(globValue, pc, ioFunctionSet).toList
    } yield WomGlobFile(path): WomFile

    files.toSet.validNel[String]
  }
}

final case class WorkflowStepInputExpression(input: WorkflowStepInput, override val cwlExpressionType: WomType, graphInputs: Set[String])(implicit parentName: ParentName) extends CwlWomExpression {

  override def sourceString = input.toString

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) = {
    (input.valueFrom, input.source) match {
      case (None, Some(WorkflowStepInputSource.String(id))) =>
        inputValues.
          get(FullyQualifiedName(id).id).
          toValidNel(s"could not find id $id in typeMap\n${inputValues.mkString("\n")}\nwhen evaluating $input.  Graph Inputs were ${graphInputs.mkString("\n")}")
      case _ => s"Could not do evaluateValue(${input.valueFrom}, ${input.source}), most likely it has not been implemented yet".invalidNel
    }
  }

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) = ???

  override def inputs = graphInputs ++ input.source.toSet.flatMap{ inputSource: InputSource => inputSource match {
    case WorkflowStepInputSource.String(s) => Set(FullyQualifiedName(s).id)
    case WorkflowStepInputSource.StringArray(sa) => sa.map(FullyQualifiedName(_).id).toSet
  }}
}


