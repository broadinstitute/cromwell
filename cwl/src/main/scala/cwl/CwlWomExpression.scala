package cwl

import cats.syntax.option._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cwl.InitialWorkDirRequirement.IwdrListingArrayEntry
import cwl.WorkflowStepInput.InputSource
import cwl.command.ParentName
import wom.expression.{IoFunctionSet, WomExpression}
import wom.types._
import wom.values._

import scala.concurrent.Await
import scala.concurrent.duration.Duration
import scala.util.Try

sealed trait CwlWomExpression extends WomExpression {

  def cwlExpressionType: WomType

  override def evaluateType(inputTypes: Map[String, WomType]): ErrorOr[WomType] = cwlExpressionType.validNel
}

case class JobPreparationExpression(expression: Expression,
                                    override val inputs: Set[String]) extends CwlWomExpression {
  val cwlExpressionType = WomAnyType

  override def sourceString = expression match {
    case Expression.ECMAScriptExpression(s) => s.value
    case Expression.ECMAScriptFunction(s) => s.value
  }

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet) = {
    val pc = ParameterContext().withInputs(inputValues, ioFunctionSet)
    expression.fold(EvaluateExpression).apply(pc).toErrorOr
  }

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) = Set.empty[WomFile].validNel
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
          Await.result(ioFunctionSet.glob(glob), Duration.Inf) match {
            case head :: Nil => WomString(head)
            case list => throw new RuntimeException(s"expecting a single File glob but instead got $list")
          }

        case _ => womValue
      }

    //CWL tells us the type this output is expected to be.  Attempt to coerce the actual output into this type.
    cwlExpressionType.coerceRawValue(globbedIfFile).toErrorOr
  }

  /*
  TODO:
   DB: It doesn't make sense to me that this function returns type WomFile but accepts a type to which it coerces.
   Wouldn't coerceTo always == WomFileType, and if not then what?
   */
  override def evaluateFiles(inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] = {

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

  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType) =
    "Programmer error: Shouldn't use WorkflowStepInputExpressions to find output files. You silly goose.".invalidNel

  override def inputs = graphInputs ++ input.source.toSet.flatMap{ inputSource: InputSource => inputSource match {
    case WorkflowStepInputSource.String(s) => Set(FullyQualifiedName(s).id)
    case WorkflowStepInputSource.StringArray(sa) => sa.map(FullyQualifiedName(_).id).toSet
  }}
}

final case class InitialWorkDirFileGeneratorExpression(entry: IwdrListingArrayEntry) extends CwlWomExpression {
  override def cwlExpressionType: WomType = WomSingleFileType
  override def sourceString: String = entry.toString

  override def evaluateValue(inputValues: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = entry match {
    case IwdrListingArrayEntry.StringDirent(content, StringOrExpression.String(entryname), _) =>
      Try(Await.result(ioFunctionSet.writeFile(entryname, content), Duration.Inf)).toErrorOr
    case _ => ??? // TODO WOM and the rest....
  }


  override def evaluateFiles(inputTypes: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
    "Programmer error: Shouldn't use InitialWorkDirRequirement listing to find output files. You silly goose.".invalidNel

  override def inputs: Set[String] = entry match {
    case IwdrListingArrayEntry.StringDirent(_, _, _) => Set.empty
    case _ => Set.empty // TODO WOM: For some cases this might need some...
  }
}
