package wom.types

import cats.instances.list._
import cats.syntax.traverse._
import cats.syntax.validated._
import common.validation.ErrorOr.ErrorOr
import common.validation.Validation._
import spray.json.JsObject
import wom.values.{WomMap, WomObject, WomObjectLike, WomValue}

case class WomCompositeType(typeMap: Map[String, WomType], structName: Option[String] = None) extends WomObjectTypeLike {

  private def validateType(values: Map[String, Any])(key: String, expectedType: WomType): ErrorOr[(String, WomValue)] = {
    (values.get(key), expectedType) match {
      case (Some(value), _) => expectedType.coerceRawValue(value).toErrorOr.map(key -> _)
      case (None, coerceTo: WomOptionalType) => (key -> coerceTo.none).validNel
      case (None, _) =>
        s"No value for field '$key' with non optional type '${expectedType.stableName}' has been provided".invalidNel
    }
  }

  override def validateAndCoerceValues(values: Map[String, Any]): ErrorOr[Map[String, WomValue]] = {
    typeMap.toList.traverse(Function.tupled(validateType(values))).map(_.toMap)
  }

  override protected def coercion = {
    case composite: WomObjectLike if isCoerceableFrom(composite.womType) => WomObject.withTypeUnsafe(composite.values, this)
    case map: WomMap if WomStringType.isCoerceableFrom(map.womType.keyType) => WomObject.withTypeUnsafe(map.value.map({ case (k, v) => k.valueString -> v }), this)
    case jsObject: JsObject => WomObject.withTypeUnsafe(jsObject.fields, this)
  }

  override def typeSpecificIsCoerceableFrom(otherType: WomType): Boolean = {
    otherType match {
      // This is as good as we can do here without the values
      case mapType: WomMapType if WomStringType.isCoerceableFrom(mapType.keyType) => true
      // Same here, it might not be coerceable but without the values we can't tell
      case WomObjectType => true
      case compositeType: WomCompositeType if (compositeType.typeMap.size == typeMap.size) &&
        compositeType.typeMap.forall({ case (key, typeValue) => typeMap.get(key).exists(_.isCoerceableFrom(typeValue)) }) =>
        true

      case _ => false
    }
  }

  override val friendlyName: String = structName.getOrElse("Object")

  override val stableName = {
    val fieldType = typeMap.map({
      case (key, value) => s"$key -> ${value.stableName}"
    }).mkString("\n")

    s"WomCompositeType {\n $fieldType \n}"
  }
}
