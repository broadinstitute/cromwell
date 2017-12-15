package cwl

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.either._
import cats.syntax.validated._
import cats.syntax.traverse._
import common.Checked
import common.validation.Checked._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import io.circe.{Json, JsonNumber, JsonObject}
import shapeless.Poly1
import wom.types._
import wom.values.{WomArray, WomBoolean, WomInteger, WomMap, WomOptionalValue, WomString, WomValue}

object CwlAnyToWomValue extends Poly1 {
  private val jsonFolder: Json.Folder[Checked[WomValue]] = new Json.Folder[Checked[WomValue]] {
    override def onNull = WomOptionalValue.none(WomNothingType).validNelCheck
    override def onBoolean(value: Boolean) = WomBoolean(value).validNelCheck
    override def onNumber(value: JsonNumber) = value.toInt.map(WomInteger.apply).toChecked(s"$value is not a valid Integer")
    override def onString(value: String) = WomString(value).validNelCheck
    // TODO POST 2.11 PLEASE map and flatMap all this :(
    override def onArray(value: Vector[Json]) = {
      value.headOption.map(_.foldWith(jsonFolder)) match {
        case Some(Right(head)) =>
          value.toList.traverse[ErrorOr, WomValue](_.foldWith(jsonFolder).toValidated) match {
            case Valid(v) => WomArray(WomArrayType(head.womType), v).validNelCheck
            case Invalid(errors) => Left(errors)
          }
        case Some(Left(errors)) => Left(errors)
        case None => WomOptionalValue.none(WomNothingType).validNelCheck
      }
    }
    override def onObject(value: JsonObject) = {
      value.values.headOption.map(_.foldWith(jsonFolder)) match {
        case Some(Right(head)) =>
          value.toList.traverse[ErrorOr, (WomValue, WomValue)]({
            case (k, v) => v.foldWith(jsonFolder) match {
              case Right(womValue) => (WomString(k), womValue).validNel
              case Left(errors) => errors.invalid
            }
          }) match {
            case Valid(v) =>  WomMap(WomMapType(WomStringType, head.womType), v.toMap).validNelCheck
            case Invalid(errors) => Left(errors)
          }
        case Some(Left(errors)) => Left(errors)
        case None => WomOptionalValue.none(WomNothingType).validNelCheck
      }
    }
  }
  
  implicit def fromFile: Case.Aux[File, WomType => Checked[WomValue]] = at[File] { file => womType =>
    // TODO POST 2.11 flatmap
    file.asWomValue.toEither match {
      case Right(value) => womType.coerceRawValue(value).toChecked
      case Left(errors) => Left(errors)
    }
  }

  implicit def fromJson: Case.Aux[Json, WomType => Checked[WomValue]] = at[Json] { json => womType =>
    json.foldWith(jsonFolder)
  }
}
