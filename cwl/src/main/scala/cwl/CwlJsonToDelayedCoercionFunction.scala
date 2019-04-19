package cwl

import cats.instances.list._
import cats.syntax.option._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import cwl.CwlCodecs._
import io.circe.{Json, JsonNumber, JsonObject}
import wom.executable.Executable.DelayedCoercionFunction
import wom.types._
import wom.values._

// With recursive types we could let circe parse it for us, but until we figure that out just parse it as Json and
// manually check for File / Directory structures
private [cwl] object CwlJsonToDelayedCoercionFunction extends Json.Folder[DelayedCoercionFunction] {
  private def simpleCoercion(value: Any)(womType: WomType): ErrorOr[WomValue] = womType.coerceRawValue(value).toErrorOr

  override def onNull: DelayedCoercionFunction = WomOptionalValue.none(_).validNel

  override def onBoolean(value: Boolean): DelayedCoercionFunction = simpleCoercion(value)

  override def onNumber(value: JsonNumber): DelayedCoercionFunction = {
    case WomFloatType => WomFloat(value.toDouble).validNel
    case WomIntegerType => value.toInt.map(WomInteger.apply).toValidNel(s"$value is not a valid Int")
    case WomLongType => value.toLong.map(WomLong.apply).toValidNel(s"$value is not a valid Long")
    case other => other.coerceRawValue(value.toString).toErrorOr
  }

  override def onString(value: String): DelayedCoercionFunction = simpleCoercion(value)

  override def onArray(value: Vector[Json]): DelayedCoercionFunction = {
    case womArrayType: WomArrayType =>
      value.toList
        .traverse(_.foldWith(this).apply(womArrayType.memberType))
        .map {
          WomArray(womArrayType, _)
        }

    case WomOptionalType(otherType) =>
      onArray(value).apply(otherType)
    case WomCoproductType(types) =>
      val attempts: List[ErrorOr[WomValue]] = types.toList.map(onArray(value)(_))
      attempts.find(_.isValid).getOrElse(attempts.sequence[ErrorOr, WomValue].map(_.head))
    case WomAnyType =>
      // Make an array of WomAny
      value.toList
        .traverse(_.foldWith(this).apply(WomAnyType))
        .map {
          WomArray(WomArrayType(WomAnyType), _)
        }
    case other => s"Cannot convert an array input value into a non array type: $other".invalidNel
  }

  override def onObject(value: JsonObject): DelayedCoercionFunction = {
    // CWL files are represented as Json objects, so this could be a file 
    case WomSingleFileType | WomMaybePopulatedFileType if value.toMap.get("class").flatMap(_.asString).contains("File") =>
      Json.fromJsonObject(value).as[File] match {
        case Left(errors) => errors.message.invalidNel
        case Right(file) => file.asWomValue
      }
    case WomMaybeListedDirectoryType | WomUnlistedDirectoryType if value.toMap.get("class").flatMap(_.asString).contains("Directory") =>
      Json.fromJsonObject(value).as[Directory] match {
        case Left(errors) => errors.message.invalidNel
        case Right(directory) => directory.asWomValue
      }
    case composite: WomCompositeType =>
      val foldedMap = value.toList.traverse[ErrorOr, (String, WomValue)]({
        case (k, v) =>
          composite.typeMap.get(k).map({ valueType =>
            v.foldWith(this).apply(valueType).map(k -> _)
          }).getOrElse(s"Input field $k is not defined in the composite input type $composite".invalidNel)
      }).map(_.toMap[String, WomValue])

      foldedMap.map(WomObject.withTypeUnsafe(_, composite))
    case WomObjectType =>
      val foldedMap = value.toList.traverse[ErrorOr, (String, WomValue)]({
        case (k, v) => v.foldWith(this).apply(WomAnyType).map(k -> _)
      }).map(_.toMap[String, WomValue])

      foldedMap.map(WomObject.apply)
    case WomOptionalType(otherType) => onObject(value).apply(otherType)
    case WomCoproductType(types) =>
      val attempts: List[ErrorOr[WomValue]] = types.toList.map(onObject(value)(_))
      //these are all Invalid, just taking head to satisfy required type of WomValue instead of List[WomValue]
      attempts.find(_.isValid).getOrElse(attempts.sequence[ErrorOr, WomValue].map(_.head))
    case other => s"Cannot convert an object value $value into a non array type: $other".invalidNel
  }
}
