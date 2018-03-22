package wdl.draft3.transforms.linking.expression.values

import cats.syntax.traverse._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import wom.format.MemorySize
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.ValueEvaluator
import wdl.model.draft3.graph.expression.ValueEvaluator.ops._
import wdl.shared.transforms.evaluation.values.EngineFunctions
import wdl4s.parser.MemoryUnit
import wom.expression.IoFunctionSet
import wom.types._
import wom.values.WomArray.WomArrayLike
import wom.values.{WomArray, WomFloat, WomInteger, WomOptionalValue, WomPair, WomSingleFile, WomString, WomValue}
import wom.types.coercion.ops._
import wom.types.coercion.defaults._
import wom.types.coercion.WomTypeCoercer

import scala.concurrent.duration.Duration
import scala.concurrent.{Await, Future}
import scala.util.{Success, Try}

object EngineFunctionEvaluators {
  implicit val readLinesFunctionEvaluator: ValueEvaluator[ReadLines] = new ValueEvaluator[ReadLines] {
    override def evaluateValue(a: ReadLines, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val readTsvFunctionEvaluator: ValueEvaluator[ReadTsv] = new ValueEvaluator[ReadTsv] {
    override def evaluateValue(a: ReadTsv, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val readMapFunctionEvaluator: ValueEvaluator[ReadMap] = new ValueEvaluator[ReadMap] {
    override def evaluateValue(a: ReadMap, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val readObjectFunctionEvaluator: ValueEvaluator[ReadObject] = new ValueEvaluator[ReadObject] {
    override def evaluateValue(a: ReadObject, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val readObjectsFunctionEvaluator: ValueEvaluator[ReadObjects] = new ValueEvaluator[ReadObjects] {
    override def evaluateValue(a: ReadObjects, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val readJsonFunctionEvaluator: ValueEvaluator[ReadJson] = new ValueEvaluator[ReadJson] {
    override def evaluateValue(a: ReadJson, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val readIntFunctionEvaluator: ValueEvaluator[ReadInt] = new ValueEvaluator[ReadInt] {
    override def evaluateValue(a: ReadInt, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val readStringFunctionEvaluator: ValueEvaluator[ReadString] = new ValueEvaluator[ReadString] {
    override def evaluateValue(a: ReadString, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val readFloatFunctionEvaluator: ValueEvaluator[ReadFloat] = new ValueEvaluator[ReadFloat] {
    override def evaluateValue(a: ReadFloat, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val readBooleanFunctionEvaluator: ValueEvaluator[ReadBoolean] = new ValueEvaluator[ReadBoolean] {
    override def evaluateValue(a: ReadBoolean, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val writeLinesFunctionEvaluator: ValueEvaluator[WriteLines] = new ValueEvaluator[WriteLines] {
    override def evaluateValue(a: WriteLines, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val writeTsvFunctionEvaluator: ValueEvaluator[WriteTsv] = new ValueEvaluator[WriteTsv] {
    override def evaluateValue(a: WriteTsv, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val writeMapFunctionEvaluator: ValueEvaluator[WriteMap] = new ValueEvaluator[WriteMap] {
    override def evaluateValue(a: WriteMap, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val writeObjectFunctionEvaluator: ValueEvaluator[WriteObject] = new ValueEvaluator[WriteObject] {
    override def evaluateValue(a: WriteObject, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val writeObjectsFunctionEvaluator: ValueEvaluator[WriteObjects] = new ValueEvaluator[WriteObjects] {
    override def evaluateValue(a: WriteObjects, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val writeJsonFunctionEvaluator: ValueEvaluator[WriteJson] = new ValueEvaluator[WriteJson] {
    override def evaluateValue(a: WriteJson, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val rangeFunctionEvaluator: ValueEvaluator[Range] = new ValueEvaluator[Range] {
    override def evaluateValue(a: Range, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      processValidatedSingleValue[WomInteger](a.param.evaluateValue(inputs, ioFunctionSet)) { integer =>
        WomArray(
          womType = WomArrayType(WomIntegerType, guaranteedNonEmpty = integer.value > 0),
          value = (0 until integer.value).map(WomInteger)
        ).validNel
      }
    }
  }

  implicit val transposeFunctionEvaluator: ValueEvaluator[Transpose] = new ValueEvaluator[Transpose] {
    override def evaluateValue(a: Transpose, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
      processValidatedSingleValue[WomArray](a.param.evaluateValue(inputs, ioFunctionSet)) { array =>
        EngineFunctions.transpose(array).toErrorOr
      }
  }

  implicit val lengthFunctionEvaluator: ValueEvaluator[Length] = new ValueEvaluator[Length] {
    override def evaluateValue(a: Length, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] =
      processValidatedSingleValue[WomArray](a.param.evaluateValue(inputs, ioFunctionSet)) { a => WomInteger(a.value.size).validNel }
  }

  implicit val flattenFunctionEvaluator: ValueEvaluator[Flatten] = new ValueEvaluator[Flatten] {
    override def evaluateValue(a: Flatten, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      def flatValues(v: WomValue): ErrorOr[Seq[WomValue]] = v match {
        case WomArrayLike(arrayLike) => arrayLike.value.validNel
        case other => s"inner item ${other.toWomString} was not an array-like".invalidNel
      }

      processValidatedSingleValue[WomArray](a.param.evaluateValue(inputs, ioFunctionSet)) { array =>
        val expandedValidation = array.value.toList.traverse[ErrorOr, Seq[WomValue]] { flatValues }
        expandedValidation map { expanded => WomArray(expanded.flatten) }
      } (coercer = WomArrayType(WomArrayType(WomAnyType)))
    }
  }

  implicit val selectFirstFunctionEvaluator: ValueEvaluator[SelectFirst] = new ValueEvaluator[SelectFirst] {
    override def evaluateValue(a: SelectFirst, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      processValidatedSingleValue[WomArray](a.param.evaluateValue(inputs, ioFunctionSet)) { array =>
        val firstValue = array.value collectFirst {
          case WomOptionalValue(_, Some(yay)) => yay
        }

        firstValue match {
          case Some(first) => first.validNel
          case None => s"select_first was called with ${array.size} empty values. We needed at least one to be filled.".invalidNel
        }
      } (coercer = WomArrayType(WomOptionalType(WomAnyType)))
    }
  }

  implicit val selectAllFunctionEvaluator: ValueEvaluator[SelectAll] = new ValueEvaluator[SelectAll] {
    override def evaluateValue(a: SelectAll, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val definedFunctionEvaluator: ValueEvaluator[Defined] = new ValueEvaluator[Defined] {
    override def evaluateValue(a: Defined, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val floorFunctionEvaluator: ValueEvaluator[Floor] = new ValueEvaluator[Floor] {
    override def evaluateValue(a: Floor, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      processValidatedSingleValue[WomFloat](a.param.evaluateValue(inputs, ioFunctionSet)) { float =>
        WomInteger(math.floor(float.value).toInt).validNel
      }
    }
  }

  implicit val ceilFunctionEvaluator: ValueEvaluator[Ceil] = new ValueEvaluator[Ceil] {
    override def evaluateValue(a: Ceil, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      processValidatedSingleValue[WomFloat](a.param.evaluateValue(inputs, ioFunctionSet)) { float =>
        WomInteger(math.ceil(float.value).toInt).validNel
      }
    }
  }

  implicit val roundFunctionEvaluator: ValueEvaluator[Round] = new ValueEvaluator[Round] {
    override def evaluateValue(a: Round, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      processValidatedSingleValue[WomFloat](a.param.evaluateValue(inputs, ioFunctionSet)) { float =>
        WomInteger(math.round(float.value).toInt).validNel
      }
    }
  }

  implicit val sizeFunctionEvaluator: ValueEvaluator[Size] = new ValueEvaluator[Size] {
    override def evaluateValue(a: Size, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      // Inner function: get the memory unit from the second (optional) parameter
      def toUnit(womValue: WomValue): ErrorOr[MemoryUnit] = Try(MemoryUnit.fromSuffix(womValue.valueString)).toErrorOr

      // Inner function: is this a file type, or an optional containing a file type?
      def isOptionalOfFileType(womType: WomType): Boolean = womType match {
        case f if WomSingleFileType.isCoerceableFrom(f) => true
        case WomOptionalType(inner) => isOptionalOfFileType(inner)
        case _ => false
      }

      // Inner function: Get the file size, allowing for unpacking of optionals and arrays
      def optionalSafeFileSize(value: WomValue): ErrorOr[Long] = value match {
        case f if f.isInstanceOf[WomSingleFile] || WomSingleFileType.isCoerceableFrom(f.womType) =>
          f.coerceToType[WomSingleFile] flatMap { file => Try(Await.result(ioFunctionSet.size(file.valueString), Duration.Inf)).toErrorOr }
        case WomOptionalValue(f, Some(o)) if isOptionalOfFileType(f) => optionalSafeFileSize(o)
        case WomOptionalValue(f, None) if isOptionalOfFileType(f) => 0l.validNel
        case WomArray(WomArrayType(womType), values) if isOptionalOfFileType(womType) => values.toList.traverse(optionalSafeFileSize).map(_.sum)
        case _ => s"The 'size' method expects a 'File', 'File?', 'Array[File]' or Array[File?] argument but instead got ${value.womType.toDisplayString}.".invalidNel
      }

      // Inner function: get the file size and convert into the requested memory unit
      def fileSize(womValue: ErrorOr[WomValue], convertTo: ErrorOr[MemoryUnit] = MemoryUnit.Bytes.validNel): ErrorOr[Double] = {
        for {
          value <- womValue
          unit <- convertTo
          fileSize <- optionalSafeFileSize(value)
        } yield MemorySize(fileSize.toDouble, MemoryUnit.Bytes).to(unit).amount
      }

      val evaluatedFileValidation: ErrorOr[WomValue] = a.file.evaluateValue(inputs, ioFunctionSet)
      a.unit map (_.evaluateValue(inputs, ioFunctionSet)) match {
        case None => fileSize(evaluatedFileValidation) map WomFloat.apply
        case Some(evaluatedUnitValidation) => fileSize(evaluatedFileValidation, evaluatedUnitValidation flatMap toUnit) map WomFloat.apply
      }
    }
  }

  implicit val basenameFunctionEvaluator: ValueEvaluator[Basename] = new ValueEvaluator[Basename] {
    override def evaluateValue(a: Basename, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      def simpleBasename(fileNameAsString: WomString) = fileNameAsString.valueString.split('/').last

      a.suffixToRemove match {
        case None => processValidatedSingleValue[WomString](a.param.evaluateValue(inputs, ioFunctionSet)) { str =>
          WomString(simpleBasename(str)).validNel
        }
        case Some(suffixToRemove) => processTwoValidatedValues[WomString, WomString](
          a.param.evaluateValue(inputs, ioFunctionSet),
          suffixToRemove.evaluateValue(inputs, ioFunctionSet)) { (name, suffix) =>
            WomString(simpleBasename(name).stripSuffix(suffix.valueString)).validNel
          }
      }
    }
  }

  implicit val zipFunctionEvaluator: ValueEvaluator[Zip] = new ValueEvaluator[Zip] {
    override def evaluateValue(a: Zip, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      processTwoValidatedValues[WomArray, WomArray](a.arg1.evaluateValue(inputs, ioFunctionSet), a.arg2.evaluateValue(inputs, ioFunctionSet)) { (arr1, arr2) =>
        if (arr1.size == arr2.size) {
          val pairs = arr1.value.zip(arr2.value) map { case (a, b) => WomPair(a, b) }
          WomArray(WomArrayType(WomPairType(arr1.arrayType.memberType, arr2.arrayType.memberType)), pairs).validNel
        } else {
          s"Mismatching array sizes for zip function: ${arr1.size} vs ${arr2.size}".invalidNel
        }
      }
    }
  }

  implicit val crossFunctionEvaluator: ValueEvaluator[Cross] = new ValueEvaluator[Cross] {
    override def evaluateValue(a: Cross, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = {
      processTwoValidatedValues[WomArray, WomArray](a.arg1.evaluateValue(inputs, ioFunctionSet), a.arg2.evaluateValue(inputs, ioFunctionSet)) { (arr1, arr2) =>
        val pairs = for {
          a <- arr1.value
          b <- arr2.value
        } yield WomPair(a, b)
        WomArray(WomArrayType(WomPairType(arr1.arrayType.memberType, arr2.arrayType.memberType)), pairs).validNel
      }
    }
  }

  implicit val prefixFunctionEvaluator: ValueEvaluator[Prefix] = new ValueEvaluator[Prefix] {
    override def evaluateValue(a: Prefix, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  implicit val subFunctionEvaluator: ValueEvaluator[Sub] = new ValueEvaluator[Sub] {
    override def evaluateValue(a: Sub, inputs: Map[String, WomValue], ioFunctionSet: IoFunctionSet): ErrorOr[WomValue] = ???
  }

  private def processValidatedSingleValue[A <: WomValue](arg: ErrorOr[WomValue])(f: A => ErrorOr[WomValue])(implicit coercer: WomTypeCoercer[A]): ErrorOr[WomValue] = {
    arg flatMap {
      case a: WomValue if a.coercionDefined[A] => a.coerceToType[A] flatMap { f.apply }
      case other => s"Expected ${coercer.toDisplayString} argument but got ${other.womType.toDisplayString}".invalidNel
    }
  }

  private def processTwoValidatedValues[A <: WomValue, B <: WomValue](arg1: ErrorOr[WomValue], arg2: ErrorOr[WomValue])
                                                                     (f: (A, B) => ErrorOr[WomValue])
                                                                     (implicit coercerA: WomTypeCoercer[A],
                                                                      coercerB: WomTypeCoercer[B]): ErrorOr[WomValue] = {
    (arg1, arg2) flatMapN {
      case (a, b) if a.coercionDefined[A] && b.coercionDefined[B] => (a.coerceToType[A], b.coerceToType[B]) flatMapN { f.apply }
      case (otherA, otherB) => s"Expected (${coercerA.toDisplayString}, ${coercerB.toDisplayString}) argument but got (${otherA.womType.toDisplayString}, ${otherB.womType.toDisplayString})".invalidNel
    }
  }
}
