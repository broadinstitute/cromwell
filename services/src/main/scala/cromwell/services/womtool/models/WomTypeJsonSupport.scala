package cromwell.services.womtool.models

import spray.json._
import wom.types._

object WomTypeJsonSupport {

  // We use `wom.types.WomType.callCachingName` instead of `wom.types.WomType.displayName` here because
  // the type hierarchy is designed for machine readability and should similarly be stable.
  implicit val womTypeEncoder: JsonWriter[WomType] = new JsonWriter[WomType] {
    final def write(a: WomType): JsValue =
      a match {
        case a: WomMapType =>
          JsObject(
            ("typeName", JsString("Map")),
            ("mapType",
             JsObject(
               ("keyType", write(a.keyType)),
               ("valueType", write(a.valueType))
             )
            )
          )
        case a: WomPairType =>
          JsObject(
            ("typeName", JsString("Pair")),
            ("pairType",
             JsObject(
               ("leftType", write(a.leftType)),
               ("rightType", write(a.rightType))
             )
            )
          )
        case a: WomArrayType =>
          JsObject(
            ("typeName", JsString("Array")),
            ("arrayType", write(a.memberType)),
            ("nonEmpty", JsBoolean(a.guaranteedNonEmpty))
          )
        case a: WomCompositeType =>
          JsObject(
            ("typeName", JsString("Object")),
            ("objectFieldTypes",
             JsArray(
               a.typeMap.toVector map { entry =>
                 JsObject(
                   ("fieldName", JsString(entry._1)),
                   ("fieldType", write(entry._2))
                 )
               }
             )
            )
          )
        case a: WomOptionalType =>
          JsObject(
            ("typeName", JsString("Optional")),
            ("optionalType", write(a.memberType))
          )
        case _ =>
          JsObject(
            ("typeName", JsString(a.stableName))
          )
      }
  }

}
