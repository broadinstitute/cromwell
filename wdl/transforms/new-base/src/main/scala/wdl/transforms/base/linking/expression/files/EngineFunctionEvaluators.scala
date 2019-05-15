package wdl.transforms.base.linking.expression.files

import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{EvaluatedValue, FileEvaluator, ValueEvaluator}
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._
import wdl.transforms.base.wdlom2wdl.WdlWriter.ops._
import wdl.transforms.base.wdlom2wdl.WdlWriterImpl.expressionElementWriter

object EngineFunctionEvaluators {

  private def evaluateToFile(forFunction: String, outermostElement: ExpressionElement, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, nestedElement: Option[ExpressionElement] = None)
                            (implicit expressionValueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {
    // `IndexAccess` `ExpressionElement`s require recursion to look for contained `IdentifierLookup` or `IdentifierMemberAccess`.
    // After all `IndexAccess` elements have been traversed the original `outermostElement` should be evaluated for files.
    val elementToExamine = nestedElement getOrElse outermostElement
    (elementToExamine match {
      case IndexAccess(expressionElement, _) => evaluateToFile(forFunction, outermostElement, inputs, ioFunctionSet, Option(expressionElement))
      // If the specified identifier is not among the inputs there are no files to delocalize from this expression.
      case IdentifierLookup(identifier) if !inputs.contains(identifier) => Set.empty[WomFile].validNel
      case IdentifierMemberAccess(first, _, _) if !inputs.contains(first) => Set.empty[WomFile].validNel
      case _ =>
        outermostElement.evaluateValue(inputs, ioFunctionSet, None) flatMap {
          case EvaluatedValue(p: WomPrimitive, _) => Set[WomFile](WomSingleFile(p.valueString)).validNel
          case other => s"Expected a primitive but got ${other.getClass.getSimpleName}".invalidNel
        }
    }).contextualizeErrors(s"predict files needed to de-localize from '${outermostElement.toWdlV1}' for $forFunction")
  }

  def singleParameterPassthroughFileEvaluator[A <: OneParamFunctionCallElement]: FileEvaluator[A] = new FileEvaluator[A] {
    override def predictFilesNeededToEvaluate(a: A,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement], valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
    }
  }

  def singleParameterEvaluateToFileFileEvaluator[A <: OneParamFunctionCallElement](functionName: String): FileEvaluator[A] = new FileEvaluator[A] {
    override def predictFilesNeededToEvaluate(a: A,
                                              inputs: Map[String, WomValue],
                                              ioFunctionSet: IoFunctionSet,
                                              coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement], valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {
      evaluateToFile(functionName, a.param, inputs, ioFunctionSet)
    }
  }

  implicit val stdoutFunctionEvaluator: FileEvaluator[StdoutElement.type] = new FileEvaluator[StdoutElement.type] {
    override def predictFilesNeededToEvaluate(a: ExpressionElement.StdoutElement.type, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val stderrFunctionEvaluator: FileEvaluator[StderrElement.type] = new FileEvaluator[StderrElement.type] {
    override def predictFilesNeededToEvaluate(a: ExpressionElement.StderrElement.type, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val readLinesFunctionEvaluator: FileEvaluator[ReadLines] = singleParameterEvaluateToFileFileEvaluator("read_lines")

  implicit val readTsvFunctionEvaluator: FileEvaluator[ReadTsv] = singleParameterEvaluateToFileFileEvaluator("read_tsv")

  implicit val readMapFunctionEvaluator: FileEvaluator[ReadMap] = singleParameterEvaluateToFileFileEvaluator("read_map")

  implicit val readObjectFunctionEvaluator: FileEvaluator[ReadObject] = singleParameterEvaluateToFileFileEvaluator("read_object")

  implicit val readObjectsFunctionEvaluator: FileEvaluator[ReadObjects] = singleParameterEvaluateToFileFileEvaluator("read_objects")

  implicit val readJsonFunctionEvaluator: FileEvaluator[ReadJson] = singleParameterEvaluateToFileFileEvaluator("read_json")

  implicit val readIntFunctionEvaluator: FileEvaluator[ReadInt] = singleParameterEvaluateToFileFileEvaluator("read_int")

  implicit val readStringFunctionEvaluator: FileEvaluator[ReadString] = singleParameterEvaluateToFileFileEvaluator("read_string")

  implicit val readFloatFunctionEvaluator: FileEvaluator[ReadFloat] = singleParameterEvaluateToFileFileEvaluator("read_float")

  implicit val readBooleanFunctionEvaluator: FileEvaluator[ReadBoolean] = singleParameterEvaluateToFileFileEvaluator("read_boolean")

  implicit val writeLinesFunctionEvaluator: FileEvaluator[WriteLines] = singleParameterPassthroughFileEvaluator

  implicit val writeTsvFunctionEvaluator: FileEvaluator[WriteTsv] = singleParameterPassthroughFileEvaluator

  implicit val writeMapFunctionEvaluator: FileEvaluator[WriteMap] = singleParameterPassthroughFileEvaluator

  implicit val writeObjectFunctionEvaluator: FileEvaluator[WriteObject] = singleParameterPassthroughFileEvaluator

  implicit val writeObjectsFunctionEvaluator: FileEvaluator[WriteObjects] = singleParameterPassthroughFileEvaluator

  implicit val writeJsonFunctionEvaluator: FileEvaluator[WriteJson] = singleParameterPassthroughFileEvaluator

  implicit val rangeFunctionEvaluator: FileEvaluator[Range] = singleParameterPassthroughFileEvaluator

  implicit val transposeFunctionEvaluator: FileEvaluator[Transpose] = singleParameterPassthroughFileEvaluator

  implicit val lengthFunctionEvaluator: FileEvaluator[Length] = singleParameterPassthroughFileEvaluator

  implicit val flattenFunctionEvaluator: FileEvaluator[Flatten] = singleParameterPassthroughFileEvaluator

  implicit val selectFirstFunctionEvaluator: FileEvaluator[SelectFirst] = singleParameterPassthroughFileEvaluator

  implicit val selectAllFunctionEvaluator: FileEvaluator[SelectAll] = singleParameterPassthroughFileEvaluator

  implicit val definedFunctionEvaluator: FileEvaluator[Defined] = singleParameterPassthroughFileEvaluator

  implicit val floorFunctionEvaluator: FileEvaluator[Floor] = singleParameterPassthroughFileEvaluator

  implicit val ceilFunctionEvaluator: FileEvaluator[Ceil] = singleParameterPassthroughFileEvaluator

  implicit val roundFunctionEvaluator: FileEvaluator[Round] = singleParameterPassthroughFileEvaluator

  implicit val globFunctionEvaluator: FileEvaluator[Glob] = new FileEvaluator[Glob] {
    override def predictFilesNeededToEvaluate(a: Glob, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] = {
      a.param.evaluateValue(inputs, ioFunctionSet, None) flatMap {
        case EvaluatedValue(p: WomPrimitive, _) => Set[WomFile](WomGlobFile(p.valueString)).validNel
        case other => s"Could not predict files to delocalize from '$a' for 'glob'. Expected a primitive but got ${other.getClass.getSimpleName}".invalidNel
      }
    }
  }

  implicit val sizeFunctionEvaluator: FileEvaluator[Size] = new FileEvaluator[Size] {
    override def predictFilesNeededToEvaluate(a: Size, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] =
      (evaluateToFile("size", a.file, inputs, ioFunctionSet): ErrorOr[Set[WomFile]],
        a.unit.fold(Set.empty[WomFile].validNel: ErrorOr[Set[WomFile]])(_.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo))) mapN { _ ++ _ }
  }

  implicit val basenameFunctionEvaluator: FileEvaluator[Basename] = new FileEvaluator[Basename] {
    override def predictFilesNeededToEvaluate(a: Basename, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val zipFunctionEvaluator: FileEvaluator[Zip] = new FileEvaluator[Zip] {
    override def predictFilesNeededToEvaluate(a: Zip, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] =
      (a.arg1.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.arg2.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ }
  }

  implicit val crossFunctionEvaluator: FileEvaluator[Cross] = new FileEvaluator[Cross] {
    override def predictFilesNeededToEvaluate(a: Cross, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] =
      (a.arg1.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.arg2.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ }
  }

  implicit val prefixFunctionEvaluator: FileEvaluator[Prefix] = new FileEvaluator[Prefix] {
    override def predictFilesNeededToEvaluate(a: Prefix, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] =
      (a.arg1.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.arg2.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ }
  }

  implicit val subFunctionEvaluator: FileEvaluator[Sub] = new FileEvaluator[Sub] {
    override def predictFilesNeededToEvaluate(a: Sub, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType)
                                             (implicit fileEvaluator: FileEvaluator[ExpressionElement],
                                              valueEvaluator: ValueEvaluator[ExpressionElement]): ErrorOr[Set[WomFile]] =
      (a.pattern.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.input.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.replace.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ ++ _ }
  }
}
