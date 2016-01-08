package wdl4s.types

import wdl4s.Call
import wdl4s.values._
import spray.json.JsObject

import scala.util.{Try, Failure, Success}

case object WdlObjectType extends WdlType {
  val toWdlString: String = "Object"

  private def handleCoercionFailures(tries: Try[_]*) = {
    val errorMessages = tries collect {
      case Failure(f) => f.getMessage
    } mkString ","

    throw new UnsupportedOperationException(s"Coercion failed: $errorMessages")
  }

  override protected def coercion = {
    case o: WdlObject => o
    case m: WdlMap if isMapCoercable(m) =>
      val coercedMap = m.value map {
        case (k, v) => toWdlString(k) -> toWdlString(v)
      } collect {
        case (Success(k), Success(v)) => k.value -> v
        case (k, v) => handleCoercionFailures(k, v)
      }

      WdlObject(coercedMap)
    case js: JsObject =>
      val coercedMap = WdlMap.coerceMap(js.fields, WdlMapType(WdlStringType, WdlAnyType)).value map {
        // get is safe because coerceMap above would have failed already if k was not coerceable to WdlString
        case (k, v) => toWdlString(k).get.value -> v
      }

      WdlObject(coercedMap)
  }

  private def toWdlString(v: WdlValue) = WdlStringType.coerceRawValue(v).map(_.asInstanceOf[WdlString])

  override def isCoerceableFrom(otherType: WdlType) = otherType match {
    case WdlObjectType => true
    case t: WdlMapType if isMapTypeCoercable(t) => true
    case _ => false
  }

  def isMapTypeCoercable(t: WdlMapType) = WdlStringType.isCoerceableFrom(t.keyType) && WdlStringType.isCoerceableFrom(t.valueType)
  def isMapCoercable(m: WdlMap) = isMapTypeCoercable(m.wdlType)
}

case class WdlCallOutputsObjectType(call: Call) extends WdlType {
  val toWdlString: String = "Object"

  override protected def coercion = {
    case o: WdlCallOutputsObject => o
  }
}
