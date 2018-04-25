package wdl.draft3.transforms.linking.expression.types

import cats.data.Validated.Valid
import cats.syntax.apply._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.draft3.elements.ExpressionElement._
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.expression.TypeEvaluator.ops._
import wom.types._
import wom.values.WomArray.WomArrayLike

object EngineFunctionEvaluators {

  implicit val readLinesFunctionEvaluator: TypeEvaluator[ReadLines] = new TypeEvaluator[ReadLines] {
    override def evaluateType(a: ReadLines, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomSingleFileType).map(_ => WomArrayType(WomStringType))
    }
  }

  implicit val readTsvFunctionEvaluator: TypeEvaluator[ReadTsv] = new TypeEvaluator[ReadTsv] {
    override def evaluateType(a: ReadTsv, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomSingleFileType).map(_ => WomArrayType(WomArrayType(WomStringType)))
    }
  }

  implicit val readMapFunctionEvaluator: TypeEvaluator[ReadMap] = new TypeEvaluator[ReadMap] {
    override def evaluateType(a: ReadMap, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomSingleFileType).map(_ => WomMapType(WomStringType, WomStringType))
    }
  }

  implicit val readObjectFunctionEvaluator: TypeEvaluator[ReadObject] = new TypeEvaluator[ReadObject] {
    override def evaluateType(a: ReadObject, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomSingleFileType).map(_ => WomObjectType)
    }
  }

  implicit val readObjectsFunctionEvaluator: TypeEvaluator[ReadObjects] = new TypeEvaluator[ReadObjects] {
    override def evaluateType(a: ReadObjects, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomSingleFileType).map(_ => WomArrayType(WomObjectType))
    }
  }

  implicit val readJsonFunctionEvaluator: TypeEvaluator[ReadJson] = new TypeEvaluator[ReadJson] {
    override def evaluateType(a: ReadJson, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomSingleFileType).map(_ => WomObjectType)
    }
  }

  implicit val readIntFunctionEvaluator: TypeEvaluator[ReadInt] = new TypeEvaluator[ReadInt] {
    override def evaluateType(a: ReadInt, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomSingleFileType).map(_ => WomIntegerType)
    }
  }

  implicit val readStringFunctionEvaluator: TypeEvaluator[ReadString] = new TypeEvaluator[ReadString] {
    override def evaluateType(a: ReadString, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomSingleFileType).map(_ => WomStringType)
    }
  }

  implicit val readFloatFunctionEvaluator: TypeEvaluator[ReadFloat] = new TypeEvaluator[ReadFloat] {
    override def evaluateType(a: ReadFloat, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomSingleFileType).map(_ => WomFloatType)
    }
  }

  implicit val readBooleanFunctionEvaluator: TypeEvaluator[ReadBoolean] = new TypeEvaluator[ReadBoolean] {
    override def evaluateType(a: ReadBoolean, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomSingleFileType).map(_ => WomBooleanType)
    }
  }

  implicit val writeLinesFunctionEvaluator: TypeEvaluator[WriteLines] = new TypeEvaluator[WriteLines] {
    override def evaluateType(a: WriteLines, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomArrayType(WomStringType)).map(_ => WomSingleFileType)
    }
  }

  implicit val writeTsvFunctionEvaluator: TypeEvaluator[WriteTsv] = new TypeEvaluator[WriteTsv] {
    override def evaluateType(a: WriteTsv, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomArrayType(WomArrayType(WomStringType))).map(_ => WomSingleFileType)
    }
  }

  implicit val writeMapFunctionEvaluator: TypeEvaluator[WriteMap] = new TypeEvaluator[WriteMap] {
    override def evaluateType(a: WriteMap, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomMapType(WomAnyType, WomAnyType)).map(_ => WomSingleFileType)
    }
  }

  implicit val writeObjectFunctionEvaluator: TypeEvaluator[WriteObject] = new TypeEvaluator[WriteObject] {
    override def evaluateType(a: WriteObject, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomObjectType).map(_ => WomSingleFileType)
    }
  }

  implicit val writeObjectsFunctionEvaluator: TypeEvaluator[WriteObjects] = new TypeEvaluator[WriteObjects] {
    override def evaluateType(a: WriteObjects, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomArrayType(WomObjectType)).map(_ => WomSingleFileType)
    }
  }

  implicit val writeJsonFunctionEvaluator: TypeEvaluator[WriteJson] = new TypeEvaluator[WriteJson] {
    override def evaluateType(a: WriteJson, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomObjectType).map(_ => WomSingleFileType)
    }
  }

  implicit val rangeFunctionEvaluator: TypeEvaluator[Range] = new TypeEvaluator[Range] {
    override def evaluateType(a: Range, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomIntegerType).map(_ => WomArrayType(WomIntegerType))
    }
  }

  implicit val transposeFunctionEvaluator: TypeEvaluator[Transpose] = new TypeEvaluator[Transpose] {
    override def evaluateType(a: Transpose, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      a.param.evaluateType(linkedValues).flatMap {
        case a @ WomArrayType(WomArrayType(_)) => a.validNel
        case foundType => s"Invalid parameter '${a.param}'. Expected 'Array[Array[_]]' but got '${foundType.toDisplayString}'".invalidNel
      }
    }
  }

  implicit val lengthFunctionEvaluator: TypeEvaluator[Length] = new TypeEvaluator[Length] {
    override def evaluateType(a: Length, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomArrayType(WomAnyType)).map(_ => WomIntegerType)
    }
  }

  implicit val flattenFunctionEvaluator: TypeEvaluator[Flatten] = new TypeEvaluator[Flatten] {
    override def evaluateType(a: Flatten, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      a.param.evaluateType(linkedValues).flatMap {
        case WomArrayType(inner @ WomArrayType(_)) => inner.validNel
        case foundType => s"Invalid parameter '${a.param}'. Expected 'Array[Array[_]]' but got '${foundType.toDisplayString}'".invalidNel
      }
    }
  }

  implicit val selectFirstFunctionEvaluator: TypeEvaluator[SelectFirst] = new TypeEvaluator[SelectFirst] {
    override def evaluateType(a: SelectFirst, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      a.param.evaluateType(linkedValues).flatMap {
        case WomArrayType(WomOptionalType(inner)) => inner.validNel
        case foundType => s"Invalid parameter '${a.param}'. Expected an array of optional values (eg 'Array[X?]') but got '${foundType.toDisplayString}'".invalidNel
      }
    }
  }

  implicit val selectAllFunctionEvaluator: TypeEvaluator[SelectAll] = new TypeEvaluator[SelectAll] {
    override def evaluateType(a: SelectAll, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      a.param.evaluateType(linkedValues).flatMap {
        case WomArrayType(WomOptionalType(inner)) => WomArrayType(inner).validNel
        case foundType => s"Invalid parameter '${a.param}'. Expected an array of optional values (eg 'Array[X?]') but got '${foundType.toDisplayString}'".invalidNel
      }
    }
  }

  implicit val definedFunctionEvaluator: TypeEvaluator[Defined] = new TypeEvaluator[Defined] {
    override def evaluateType(a: Defined, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomOptionalType(WomAnyType)).map(_ => WomBooleanType)
    }
  }

  implicit val floorFunctionEvaluator: TypeEvaluator[Floor] = new TypeEvaluator[Floor] {
    override def evaluateType(a: Floor, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomFloatType).map(_ => WomIntegerType)
    }
  }

  implicit val ceilFunctionEvaluator: TypeEvaluator[Ceil] = new TypeEvaluator[Ceil] {
    override def evaluateType(a: Ceil, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomFloatType).map(_ => WomIntegerType)
    }
  }

  implicit val roundFunctionEvaluator: TypeEvaluator[Round] = new TypeEvaluator[Round] {
    override def evaluateType(a: Round, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      validateParamType(a.param, linkedValues, WomFloatType).map(_ => WomIntegerType)
    }
  }

  implicit val globFunctionTypeEvaluator: TypeEvaluator[Glob] = (a, linkedValues) => {
    validateParamType(a.param, linkedValues, WomStringType).map(_ => WomArrayType(WomSingleFileType))
  }

  implicit val sizeFunctionEvaluator: TypeEvaluator[Size] = new TypeEvaluator[Size] {
    override def evaluateType(a: Size, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      val validatedSecondArg: ErrorOr[Unit] = a.secondParam match {
        case None => ().validNel
        case Some(arg) => validateParamType(arg, linkedValues, WomStringType)
      }
      (validateParamType(a.firstParam, linkedValues, WomSingleFileType),
        validatedSecondArg) mapN { (_, _) => WomFloatType }
    }
  }

  implicit val basenameFunctionEvaluator: TypeEvaluator[Basename] = new TypeEvaluator[Basename] {
    override def evaluateType(a: Basename, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      val validatedSecondArg: ErrorOr[Unit] = a.secondParam match {
        case None => ().validNel
        case Some(arg) => validateParamType(arg, linkedValues, WomStringType)
      }
      (validateParamType(a.firstParam, linkedValues, WomSingleFileType),
        validatedSecondArg) mapN { (_, _) => WomStringType }
    }
  }

  implicit val zipFunctionEvaluator: TypeEvaluator[Zip] = new TypeEvaluator[Zip] {
    override def evaluateType(a: Zip, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      (a.arg1.evaluateType(linkedValues), a.arg2.evaluateType(linkedValues)) match {
        case (Valid(WomArrayLike(left)), Valid(WomArrayLike(right))) => WomPairType(left.arrayType.memberType, right.arrayType.memberType).validNel
        case (Valid(otherLeft), Valid(WomArrayLike(_))) => s"Invalid left parameter '${a.arg1}'. Expected Array type but got '${otherLeft.toDisplayString}'".invalidNel
        case (Valid(WomArrayLike(_)), Valid(otherRight)) => s"Invalid right parameter '${a.arg2}'. Expected Array type but got '${otherRight.toDisplayString}'".invalidNel
        // One or more are invalid, so mapN function won't actually ever run:
        case (otherLeft, otherRight) => (otherLeft, otherRight) mapN { (_, _) => WomNothingType }
      }
    }
  }

  implicit val crossFunctionEvaluator: TypeEvaluator[Cross] = new TypeEvaluator[Cross] {
    override def evaluateType(a: Cross, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      (a.arg1.evaluateType(linkedValues), a.arg2.evaluateType(linkedValues)) match {
        case (Valid(WomArrayLike(left)), Valid(WomArrayLike(right))) => WomPairType(left.arrayType.memberType, right.arrayType.memberType).validNel
        case (Valid(otherLeft), Valid(WomArrayLike(_))) => s"Invalid left parameter '${a.arg1}'. Expected Array type but got '${otherLeft.toDisplayString}'".invalidNel
        case (Valid(WomArrayLike(_)), Valid(otherRight)) => s"Invalid right parameter '${a.arg2}'. Expected Array type but got '${otherRight.toDisplayString}'".invalidNel
        // One or more are invalid, so mapN function won't actually ever run:
        case (otherLeft, otherRight) => (otherLeft, otherRight) mapN { (_, _) => WomNothingType }
      }
    }
  }

  implicit val prefixFunctionEvaluator: TypeEvaluator[Prefix] = new TypeEvaluator[Prefix] {
    override def evaluateType(a: Prefix, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {
      (validateParamType(a.prefix, linkedValues, WomStringType),
        validateParamType(a.array, linkedValues, WomArrayType(WomStringType))
      ) mapN { (_, _) => WomArrayType(WomStringType) }
    }
  }

  implicit val subFunctionEvaluator: TypeEvaluator[Sub] = new TypeEvaluator[Sub] {
    override def evaluateType(a: Sub, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle]): ErrorOr[WomType] = {

      (validateParamType(a.input, linkedValues, WomSingleFileType),
        validateParamType(a.pattern, linkedValues, WomSingleFileType),
        validateParamType(a.replace, linkedValues, WomSingleFileType)) mapN { (_, _, _) => WomStringType }
    }
  }

  private def validateParamType(param: ExpressionElement, linkedValues: Map[UnlinkedConsumedValueHook, GeneratedValueHandle], expectedType: WomType): ErrorOr[Unit] = {
    param.evaluateType(linkedValues).flatMap { foundType =>
      if (expectedType.isCoerceableFrom(foundType)) { ().validNel } else { s"Invalid parameter '$param'. Expected '${expectedType.toDisplayString}' but got '${foundType.toDisplayString}'".invalidNel }
    }
  }
}
