package wdl.draft3.transforms.linking.expression.files

import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import wdl.model.draft3.elements.ExpressionElement
import wdl.draft3.transforms.linking.expression.values._
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.{EvaluatedValue, FileEvaluator}
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.model.draft3.graph.expression.FileEvaluator.ops._
import wom.expression.IoFunctionSet
import wom.types._
import wom.values._
import wdl.draft3.transforms.wdlom2wdl.WdlWriter.ops._
import wdl.draft3.transforms.wdlom2wdl.WdlWriterImpl.expressionElementWriter

object EngineFunctionEvaluators {

  private def evaluateToFile(forFunction: String, a: ExpressionElement, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[Set[WomFile]] = {
    (a match {
      case _: IdentifierLookup | _: IdentifierMemberAccess | _: IndexAccess => Set.empty[WomFile].validNel
      case _ =>
        a.evaluateValue(inputs, ioFunctionSet, None) flatMap {
          case EvaluatedValue(p: WomPrimitive, _) => Set[WomFile](WomSingleFile(p.valueString)).validNel
          case other => s"Expected a primitive but got ${other.getClass.getSimpleName}".invalidNel
        }
    }).contextualizeErrors(s"predict files needed to de-localize from '${a.toWdlV1}' for $forFunction")
  }

  implicit val stdoutFunctionEvaluator: FileEvaluator[StdoutElement.type] = new FileEvaluator[StdoutElement.type] {
    override def predictFilesNeededToEvaluate(a: ExpressionElement.StdoutElement.type, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val stderrFunctionEvaluator: FileEvaluator[StderrElement.type] = new FileEvaluator[StderrElement.type] {
    override def predictFilesNeededToEvaluate(a: ExpressionElement.StderrElement.type, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val readLinesFunctionEvaluator: FileEvaluator[ReadLines] = new FileEvaluator[ReadLines] {
    override def predictFilesNeededToEvaluate(a: ReadLines, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      evaluateToFile("read_lines", a.param, inputs, ioFunctionSet)
  }

  implicit val readTsvFunctionEvaluator: FileEvaluator[ReadTsv] = new FileEvaluator[ReadTsv] {
    override def predictFilesNeededToEvaluate(a: ReadTsv, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      evaluateToFile("read_tsv", a.param, inputs, ioFunctionSet)
  }

  implicit val readMapFunctionEvaluator: FileEvaluator[ReadMap] = new FileEvaluator[ReadMap] {
    override def predictFilesNeededToEvaluate(a: ReadMap, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      evaluateToFile("read_map", a.param, inputs, ioFunctionSet)
  }

  implicit val readObjectFunctionEvaluator: FileEvaluator[ReadObject] = new FileEvaluator[ReadObject] {
    override def predictFilesNeededToEvaluate(a: ReadObject, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      evaluateToFile("read_object", a.param, inputs, ioFunctionSet)
  }

  implicit val readObjectsFunctionEvaluator: FileEvaluator[ReadObjects] = new FileEvaluator[ReadObjects] {
    override def predictFilesNeededToEvaluate(a: ReadObjects, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      evaluateToFile("read_objects", a.param, inputs, ioFunctionSet)
  }

  implicit val readJsonFunctionEvaluator: FileEvaluator[ReadJson] = new FileEvaluator[ReadJson] {
    override def predictFilesNeededToEvaluate(a: ReadJson, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      evaluateToFile("read_json", a.param, inputs, ioFunctionSet)
  }

  implicit val readIntFunctionEvaluator: FileEvaluator[ReadInt] = new FileEvaluator[ReadInt] {
    override def predictFilesNeededToEvaluate(a: ReadInt, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      evaluateToFile("read_int", a.param, inputs, ioFunctionSet)
  }

  implicit val readStringFunctionEvaluator: FileEvaluator[ReadString] = new FileEvaluator[ReadString] {
    override def predictFilesNeededToEvaluate(a: ReadString, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      evaluateToFile("read_string", a.param, inputs, ioFunctionSet)
  }

  implicit val readFloatFunctionEvaluator: FileEvaluator[ReadFloat] = new FileEvaluator[ReadFloat] {
    override def predictFilesNeededToEvaluate(a: ReadFloat, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      evaluateToFile("read_float", a.param, inputs, ioFunctionSet)
  }

  implicit val readBooleanFunctionEvaluator: FileEvaluator[ReadBoolean] = new FileEvaluator[ReadBoolean] {
    override def predictFilesNeededToEvaluate(a: ReadBoolean, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      evaluateToFile("read_boolean", a.param, inputs, ioFunctionSet)
  }

  implicit val writeLinesFunctionEvaluator: FileEvaluator[WriteLines] = new FileEvaluator[WriteLines] {
    override def predictFilesNeededToEvaluate(a: WriteLines, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val writeTsvFunctionEvaluator: FileEvaluator[WriteTsv] = new FileEvaluator[WriteTsv] {
    override def predictFilesNeededToEvaluate(a: WriteTsv, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val writeMapFunctionEvaluator: FileEvaluator[WriteMap] = new FileEvaluator[WriteMap] {
    override def predictFilesNeededToEvaluate(a: WriteMap, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val writeObjectFunctionEvaluator: FileEvaluator[WriteObject] = new FileEvaluator[WriteObject] {
    override def predictFilesNeededToEvaluate(a: WriteObject, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val writeObjectsFunctionEvaluator: FileEvaluator[WriteObjects] = new FileEvaluator[WriteObjects] {
    override def predictFilesNeededToEvaluate(a: WriteObjects, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val writeJsonFunctionEvaluator: FileEvaluator[WriteJson] = new FileEvaluator[WriteJson] {
    override def predictFilesNeededToEvaluate(a: WriteJson, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val rangeFunctionEvaluator: FileEvaluator[Range] = new FileEvaluator[Range] {
    override def predictFilesNeededToEvaluate(a: Range, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }

  implicit val transposeFunctionEvaluator: FileEvaluator[Transpose] = new FileEvaluator[Transpose] {
    override def predictFilesNeededToEvaluate(a: Transpose, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }

  implicit val lengthFunctionEvaluator: FileEvaluator[Length] = new FileEvaluator[Length] {
    override def predictFilesNeededToEvaluate(a: Length, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }

  implicit val flattenFunctionEvaluator: FileEvaluator[Flatten] = new FileEvaluator[Flatten] {
    override def predictFilesNeededToEvaluate(a: Flatten, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }

  implicit val selectFirstFunctionEvaluator: FileEvaluator[SelectFirst] = new FileEvaluator[SelectFirst] {
    override def predictFilesNeededToEvaluate(a: SelectFirst, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }

  implicit val selectAllFunctionEvaluator: FileEvaluator[SelectAll] = new FileEvaluator[SelectAll] {
    override def predictFilesNeededToEvaluate(a: SelectAll, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }

  implicit val definedFunctionEvaluator: FileEvaluator[Defined] = new FileEvaluator[Defined] {
    override def predictFilesNeededToEvaluate(a: Defined, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }

  implicit val floorFunctionEvaluator: FileEvaluator[Floor] = new FileEvaluator[Floor] {
    override def predictFilesNeededToEvaluate(a: Floor, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }

  implicit val ceilFunctionEvaluator: FileEvaluator[Ceil] = new FileEvaluator[Ceil] {
    override def predictFilesNeededToEvaluate(a: Ceil, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }

  implicit val roundFunctionEvaluator: FileEvaluator[Round] = new FileEvaluator[Round] {
    override def predictFilesNeededToEvaluate(a: Round, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      a.param.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)
  }

  implicit val globFunctionEvaluator: FileEvaluator[Glob] = (a, inputs, ioFunctionSet, _) => {
    a.param.evaluateValue(inputs, ioFunctionSet, None) flatMap {
      case EvaluatedValue(p: WomPrimitive, _) => Set[WomFile](WomGlobFile(p.valueString)).validNel
      case other => s"Could not predict files to delocalize from '$a' for 'glob'. Expected a primitive but got ${other.getClass.getSimpleName}".invalidNel
    }
  }

  implicit val sizeFunctionEvaluator: FileEvaluator[Size] = new FileEvaluator[Size] {
    override def predictFilesNeededToEvaluate(a: Size, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      (evaluateToFile("size", a.file, inputs, ioFunctionSet): ErrorOr[Set[WomFile]],
        a.unit.fold(Set.empty[WomFile].validNel: ErrorOr[Set[WomFile]])(_.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo))) mapN { _ ++ _ }
  }

  implicit val basenameFunctionEvaluator: FileEvaluator[Basename] = new FileEvaluator[Basename] {
    override def predictFilesNeededToEvaluate(a: Basename, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      Set.empty[WomFile].validNel
  }

  implicit val zipFunctionEvaluator: FileEvaluator[Zip] = new FileEvaluator[Zip] {
    override def predictFilesNeededToEvaluate(a: Zip, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      (a.arg1.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.arg2.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ }
  }

  implicit val crossFunctionEvaluator: FileEvaluator[Cross] = new FileEvaluator[Cross] {
    override def predictFilesNeededToEvaluate(a: Cross, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      (a.arg1.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.arg2.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ }
  }

  implicit val prefixFunctionEvaluator: FileEvaluator[Prefix] = new FileEvaluator[Prefix] {
    override def predictFilesNeededToEvaluate(a: Prefix, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      (a.arg1.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.arg2.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ }
  }

  implicit val subFunctionEvaluator: FileEvaluator[Sub] = new FileEvaluator[Sub] {
    override def predictFilesNeededToEvaluate(a: Sub, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet, coerceTo: WomType): ErrorOr[Set[WomFile]] =
      (a.pattern.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.input.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo),
        a.replace.evaluateFilesNeededToEvaluate(inputs, ioFunctionSet, coerceTo)) mapN { _ ++ _ ++ _ }
  }
}
