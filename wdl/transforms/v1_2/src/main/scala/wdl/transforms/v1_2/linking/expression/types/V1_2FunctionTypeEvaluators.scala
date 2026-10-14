package wdl.transforms.v1_2.linking.expression.types

import cats.syntax.apply._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr._
import wdl.model.draft3.elements.ExpressionElement
import wdl.model.draft3.elements.ExpressionElement.{
  Equals,
  Length,
  NotEquals,
  PrimitiveLiteralExpressionElement,
  Size,
  WriteTsv
}
import wdl.model.draft3.graph.expression.TypeEvaluator
import wdl.model.draft3.graph.{GeneratedValueHandle, UnlinkedConsumedValueHook}
import wdl.model.v1_2.elements.V1_2ExpressionElement._
import wom.types._
import wom.values.WomBoolean

/**
 * Type evaluation for the operators and functions added in WDL 1.2, and for the WDL 1.0/1.1 functions whose accepted
 * argument types were extended in WDL 1.2.
 */
object V1_2FunctionTypeEvaluators {

  type LinkedValues = Map[UnlinkedConsumedValueHook, GeneratedValueHandle]

  def evaluateType(a: V1_2FunctionCallElement, linkedValues: LinkedValues, typeAliases: Map[String, WomType])(implicit
    expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
  ): ErrorOr[WomType] = {
    def typeOf(e: ExpressionElement): ErrorOr[WomType] =
      expressionTypeEvaluator.evaluateType(e, linkedValues, typeAliases)

    def expect(e: ExpressionElement, expected: WomType, functionName: String): ErrorOr[WomType] = typeOf(e) flatMap {
      case found if expected.isCoerceableFrom(found) => found.validNel
      case found =>
        s"Invalid argument to '$functionName'. Expected '${expected.stableName}' but got '${found.stableName}'".invalidNel
    }

    a match {
      case Exponentiation(left, right) =>
        (typeOf(left), typeOf(right)) flatMapN {
          case (WomIntegerType, WomIntegerType) => WomIntegerType.validNel
          case (l, r) if isNumeric(l) && isNumeric(r) => WomFloatType.validNel
          case (l, r) =>
            s"Cannot evaluate '${l.stableName} ** ${r.stableName}'. Both operands must be numeric".invalidNel
        }

      case Find(input, pattern) =>
        (expect(input, WomStringType, "find"), expect(pattern, WomStringType, "find")) mapN { (_, _) =>
          WomOptionalType(WomStringType)
        }
      case Matches(input, pattern) =>
        (expect(input, WomStringType, "matches"), expect(pattern, WomStringType, "matches")) mapN { (_, _) =>
          WomBooleanType
        }

      case JoinPaths(paths, None) => expect(paths, WomArrayType(WomStringType), "join_paths").map(_ => WomStringType)
      case JoinPaths(base, Some(relative)) =>
        val relativeType: ErrorOr[WomType] = typeOf(relative) flatMap {
          case t if WomStringType.isCoerceableFrom(t) || WomArrayType(WomStringType).isCoerceableFrom(t) => t.validNel
          case t =>
            s"Invalid second argument to 'join_paths'. Expected 'String' or 'Array[String]+' but got '${t.stableName}'".invalidNel
        }
        val baseType: ErrorOr[WomType] = typeOf(base) flatMap {
          case t if WomUnlistedDirectoryType.isCoerceableFrom(t) || WomStringType.isCoerceableFrom(t) => t.validNel
          case t => s"Invalid first argument to 'join_paths'. Expected 'Directory' but got '${t.stableName}'".invalidNel
        }
        (baseType, relativeType) mapN { (_, _) => WomStringType }

      case ReadTsvWithHeader(file, header, fieldNames) =>
        (expect(file, WomSingleFileType, "read_tsv"),
         expect(header, WomBooleanType, "read_tsv"),
         fieldNames.traverse(expect(_, WomArrayType(WomStringType), "read_tsv"))
        ) mapN { (_, _, names) =>
          header match {
            // `read_tsv(file, false)` is the same as `read_tsv(file)`
            case PrimitiveLiteralExpressionElement(WomBoolean(false)) if names.isEmpty =>
              WomArrayType(WomArrayType(WomStringType))
            case _ => WomArrayType(WomObjectType)
          }
        }

      case Contains(array, value) =>
        (typeOf(array), typeOf(value)) flatMapN {
          case (WomArrayType(memberType), valueType) if memberType.isCoerceableFrom(valueType) =>
            WomBooleanType.validNel
          case (WomArrayType(memberType), valueType) =>
            s"Cannot call 'contains' with a '${valueType.stableName}' value on an array of '${memberType.stableName}'".invalidNel
          case (other, _) =>
            s"Invalid first argument to 'contains'. Expected an array but got '${other.stableName}'".invalidNel
        }

      case Chunk(array, size) =>
        (typeOf(array), expect(size, WomIntegerType, "chunk")) flatMapN {
          case (arrayType: WomArrayType, _) => WomArrayType(WomArrayType(arrayType.memberType)).validNel
          case (other, _) =>
            s"Invalid first argument to 'chunk'. Expected an array but got '${other.stableName}'".invalidNel
        }

      case SelectFirstWithDefault(array, default) =>
        (typeOf(array), typeOf(default)) flatMapN {
          // The array may be empty, e.g. `select_first([], 5)`
          case (WomArrayType(WomNothingType | WomOptionalType(WomNothingType)), defaultType) => defaultType.validNel
          // The array members have no common type more specific than Any, e.g. `select_first([1, None], 5)`
          case (WomArrayType(WomAnyType | WomOptionalType(WomAnyType)), defaultType) => defaultType.validNel
          case (WomArrayType(WomOptionalType(inner)), defaultType) if inner.isCoerceableFrom(defaultType) =>
            inner.validNel
          case (WomArrayType(inner), defaultType) if inner.isCoerceableFrom(defaultType) => inner.validNel
          case (WomArrayType(inner), defaultType) =>
            s"The default value of 'select_first' has type '${defaultType.stableName}' but the array contains '${inner.stableName}'".invalidNel
          case (other, _) =>
            s"Invalid first argument to 'select_first'. Expected an array of optional values but got '${other.stableName}'".invalidNel
        }

      case ContainsKey(collection, key) =>
        (typeOf(collection), typeOf(key)) flatMapN {
          case (WomMapType(mapKeyType, _), keyType) if mapKeyType.isCoerceableFrom(keyType) => WomBooleanType.validNel
          case (WomMapType(_, _) | _: WomCompositeType | WomObjectType, keyType)
              if WomArrayType(WomStringType).isCoerceableFrom(keyType) =>
            WomBooleanType.validNel
          case (_: WomCompositeType | WomObjectType, keyType) if WomStringType.isCoerceableFrom(keyType) =>
            WomBooleanType.validNel
          case (collectionType, keyType) =>
            s"Cannot call 'contains_key' on '${collectionType.stableName}' with a key of type '${keyType.stableName}'".invalidNel
        }

      case Values(map) =>
        typeOf(map) flatMap {
          case WomMapType(_, valueType) => WomArrayType(valueType).validNel
          case other => s"Cannot invoke 'values' on type '${other.stableName}'. Expected a map".invalidNel
        }

      case WriteTsvWithHeader(array, header, columnNames) =>
        (validateWriteTsvArray(typeOf(array)),
         expect(header, WomBooleanType, "write_tsv"),
         columnNames.traverse(expect(_, WomArrayType(WomStringType), "write_tsv"))
        ) mapN { (_, _, _) => WomSingleFileType }

      case EnvironmentVariableValue(value) => typeOf(value).map(_ => WomStringType)
    }
  }

  // WDL 1.2: `write_tsv` also accepts an array of structs
  def evaluateWriteTsvType(a: WriteTsv, linkedValues: LinkedValues, typeAliases: Map[String, WomType])(implicit
    expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
  ): ErrorOr[WomType] =
    validateWriteTsvArray(expressionTypeEvaluator.evaluateType(a.param, linkedValues, typeAliases))
      .map(_ => WomSingleFileType)

  private def validateWriteTsvArray(arrayType: ErrorOr[WomType]): ErrorOr[Unit] = arrayType flatMap {
    case WomArrayType(_: WomCompositeType) => ().validNel
    case t if WomArrayType(WomArrayType(WomStringType)).isCoerceableFrom(t) => ().validNel
    case t =>
      s"Invalid argument to 'write_tsv'. Expected 'Array[Array[String]]' or an array of structs but got '${t.stableName}'".invalidNel
  }

  /**
   * WDL 1.2: either or both operands of `==` and `!=` may be optional, and `None` is equal only to itself.
   * Returns `None` when neither operand is optional, so that the WDL 1.0/1.1 rules apply.
   */
  def evaluateOptionalEqualityType(left: ExpressionElement,
                                   right: ExpressionElement,
                                   linkedValues: LinkedValues,
                                   typeAliases: Map[String, WomType]
  )(implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): Option[ErrorOr[WomType]] = {
    val types = (expressionTypeEvaluator.evaluateType(left, linkedValues, typeAliases),
                 expressionTypeEvaluator.evaluateType(right, linkedValues, typeAliases)
    )
    types match {
      case (cats.data.Validated.Valid(l), cats.data.Validated.Valid(r)) if isOptional(l) || isOptional(r) =>
        val (innerLeft, innerRight) = (unpackOptional(l), unpackOptional(r))
        if (
          innerLeft == WomNothingType || innerRight == WomNothingType ||
          innerLeft.isCoerceableFrom(innerRight) || innerRight.isCoerceableFrom(innerLeft)
        ) Option(WomBooleanType.validNel)
        else Option(s"Cannot compare '${l.stableName}' with '${r.stableName}'".invalidNel)
      case _ => None
    }
  }

  def evaluateEqualsType(a: Equals, linkedValues: LinkedValues, typeAliases: Map[String, WomType])(
    fallback: => ErrorOr[WomType]
  )(implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] =
    evaluateOptionalEqualityType(a.left, a.right, linkedValues, typeAliases).getOrElse(fallback)

  def evaluateNotEqualsType(a: NotEquals, linkedValues: LinkedValues, typeAliases: Map[String, WomType])(
    fallback: => ErrorOr[WomType]
  )(implicit expressionTypeEvaluator: TypeEvaluator[ExpressionElement]): ErrorOr[WomType] =
    evaluateOptionalEqualityType(a.left, a.right, linkedValues, typeAliases).getOrElse(fallback)

  private def isOptional(womType: WomType): Boolean = womType match {
    case _: WomOptionalType | WomNothingType => true
    case _ => false
  }

  private def unpackOptional(womType: WomType): WomType = womType match {
    case WomOptionalType(inner) => unpackOptional(inner)
    case other => other
  }

  // WDL 1.2: `length` also accepts a Map, an Object or a String
  def evaluateLengthType(a: Length, linkedValues: LinkedValues, typeAliases: Map[String, WomType])(implicit
    expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
  ): ErrorOr[WomType] =
    expressionTypeEvaluator.evaluateType(a.param, linkedValues, typeAliases) flatMap {
      case _: WomArrayType | _: WomMapType | WomObjectType | WomStringType => WomIntegerType.validNel
      case other =>
        s"Cannot invoke 'length' on type '${other.stableName}'. Expected an array, map, object or string".invalidNel
    }

  // WDL 1.2: `size` also accepts a Directory, or any compound value containing files or directories
  def evaluateSizeType(a: Size, linkedValues: LinkedValues, typeAliases: Map[String, WomType])(implicit
    expressionTypeEvaluator: TypeEvaluator[ExpressionElement]
  ): ErrorOr[WomType] = {
    val unitValidation: ErrorOr[Unit] =
      a.secondParam.traverse(p => expressionTypeEvaluator.evaluateType(p, linkedValues, typeAliases)) flatMap {
        case Some(t) if !WomStringType.isCoerceableFrom(t) =>
          s"Invalid unit argument to 'size'. Expected 'String' but got '${t.stableName}'".invalidNel
        case _ => ().validNel
      }
    val fileValidation: ErrorOr[Unit] =
      expressionTypeEvaluator.evaluateType(a.firstParam, linkedValues, typeAliases) flatMap {
        case t if mayContainFiles(t) => ().validNel
        case t =>
          s"Cannot invoke 'size' on type '${t.stableName}'. Expected a File, Directory, or a compound value containing them".invalidNel
      }
    (unitValidation, fileValidation) mapN { (_, _) => WomFloatType }
  }

  private def isNumeric(womType: WomType): Boolean = womType match {
    case WomIntegerType | WomFloatType => true
    case WomOptionalType(inner) => isNumeric(inner)
    case _ => false
  }

  private def mayContainFiles(womType: WomType): Boolean = womType match {
    case WomSingleFileType | WomUnlistedDirectoryType | WomStringType | WomNothingType => true
    case WomOptionalType(inner) => mayContainFiles(inner)
    case WomArrayType(inner) => mayContainFiles(inner)
    case WomMapType(keyType, valueType) => mayContainFiles(keyType) || mayContainFiles(valueType)
    case WomPairType(left, right) => mayContainFiles(left) || mayContainFiles(right)
    case WomCompositeType(members, _) => members.values.exists(mayContainFiles)
    case WomObjectType => true
    case _ => false
  }
}
