package wom.types.coercion

import cats.syntax.validated._
import common.validation.Validation._
import common.validation.ErrorOr.ErrorOr
import common.validation.ErrorOr._
import wom.types._
import wom.values._

import scala.language.implicitConversions
import scala.reflect.ClassTag

package object defaults {
  implicit val womIntegerCoercer: WomTypeCoercer[WomInteger] = defaultCoercionForType[WomInteger](WomIntegerType)
  implicit val womFloatCoercer: WomTypeCoercer[WomFloat] = defaultCoercionForType[WomFloat](WomFloatType)
  implicit val womStringCoercer: WomTypeCoercer[WomString] = defaultCoercionForType[WomString](WomStringType)
  implicit val womSingleFileCoercer: WomTypeCoercer[WomSingleFile] = defaultCoercionForType[WomSingleFile](WomSingleFileType)

  implicit val womArrayOfAnyCoercer = defaultCoercionForType[WomArray](WomArrayType(WomAnyType))
  implicit def womArrayTypeCoercer(arrayType: WomArrayType): WomTypeCoercer[WomArray] = defaultCoercionForType[WomArray](arrayType)

  private def defaultCoercionForType[A](typeObject: WomType)(implicit classTag: ClassTag[A]): WomTypeCoercer[A] = new WomTypeCoercer[A] {
    override def coerceToType(any: Any): ErrorOr[A] = typeObject.coerceRawValue(any).toErrorOr flatMap {
      case womValue: A => womValue.validNel
      case other => s"Bad coercion in ${getClass.getSimpleName}! Coercion should have created ${typeObject.toDisplayString} but instead created ${other.womType.toDisplayString}".invalidNel
    }
    override def coercionDefined(any: Any): Boolean = typeObject.coercionDefined(any)
    override def toDisplayString: String = typeObject.toDisplayString
  }
}
