package wom.util


import common.validation.ErrorOr._
import org.mozilla.javascript.{Context, ScriptableObject}
import wom.values._
import mouse.all._

object JsUtil {
  val encoder = new JsEncoder
  def writeValue(value: Js)(context: Context, scope: ScriptableObject): AnyRef =
    value match {
      case JsObject(fields) => {
        val newObj = context.newObject(scope)

        fields.toList.foreach{
          case (name, value) =>
            val newerObj = writeValue(value)(context, scope)
            ScriptableObject.putProperty(newObj, name, newerObj)
        }
        newObj
      }

      case JsArray(array) =>
        val newObj = context.newArray(scope, array.length)

        array.toList.zipWithIndex.foreach {
          case (js, index) =>

            val newerObj = writeValue(js)(context, scope)
            ScriptableObject.putProperty(newObj, index, newerObj)

        }
        newObj

      case JsField(obj) =>
        obj
    }

  def evalRaw[A](expr: String)(block: (Context, ScriptableObject) => A): AnyRef = {
    val context = Context.enter()
    try {
      context.setLanguageVersion(Context.VERSION_1_8)
      val scope = context.initStandardObjects
      block(context, scope)
      context.evaluateString(scope, expr, "not sure where this is used yet", 1, null)
    } finally {
      Context.exit()
    }
  }

  sealed trait Js

  case class JsObject(fields: Map[String, Js]) extends Js
  case class JsArray(array: Array[Js]) extends Js
  case class JsField(anyRef: AnyRef) extends Js


  /**
    * Evaluates a javascript expression using maps of WOM values.
    *
    * TODO: Once custom types are supported as WomValue, this custom method won't be required.
    *
    * @param expr      The javascript expression.
    * @param rawValues A map filled with WOM values.
    * @param mapValues A map of maps filled with WOM values of various types.
    * @param encoder   Encodes wom values to javascript.
    * @param decoder   Decodes wom values from javascript.
    * @return The result of the expression.
    */
  def evalStructish(expr: String,
                    rawValues: (String, WomValue),
                    mapValues: Map[String, Map[String, WomValue]] = Map.empty,
                    encoder: JsEncoder = new JsEncoder,
                    decoder: JsDecoder = new JsDecoder): ErrorOr[WomValue] = {
    evalRaw(expr) { (context, scope) =>

      val (key, value) = rawValues

      val jsValue = encoder.encode(value)
      val field = writeValue(jsValue)(context, scope)
      ScriptableObject.putProperty(scope, key, field)


      val jsMap = mapValues.mapValues{ _.mapValues(encoder.encode) }

      jsMap.map {
        case (scopeId, nestedMap) =>

          val newObj = context.newObject(scope)
          nestedMap.toList.map{
            case (key, value) =>
              val newerObj = writeValue(value)(context, scope)
              ScriptableObject.putProperty(newObj, key, newerObj)
          }
          ScriptableObject.putProperty(scope, scopeId, newObj)
      }

    } |> decoder.decode
  }
}
