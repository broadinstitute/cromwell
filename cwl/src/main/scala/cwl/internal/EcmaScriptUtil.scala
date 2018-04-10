package cwl.internal

import common.validation.ErrorOr._
import common.validation.Validation._
import org.mozilla.javascript._
import wom.values._

import scala.concurrent.duration._
import scala.util.Try

/**
  * This implementation depends on Mozilla Rhino.
  *
  * Previously we attempted to use Nashorn which is built into the JDK, but we encountered at least 2 situations where
  * it didn't work and we found no workarounds to satisfy the use cases.  Namely, JSON.stringify of a Map and calling "sort" on an array.
  */
//noinspection VariablePatternShadow
object EcmaScriptUtil {
  def writeValue(value: ECMAScriptVariable)(context: Context, scope: Scriptable): AnyRef =
    value match {
      case ESObject(fields) =>
        val newObj = context.newObject(scope)

        fields.toList.foreach{
          case (name, value) =>
            val newerObj = writeValue(value)(context, scope)
            ScriptableObject.putProperty(newObj, name, newerObj)
        }
        newObj

      case ESArray(array) =>
        val newObj = context.newArray(scope, array.length)

        array.toList.zipWithIndex.foreach {
          case (js, index) =>

            val newerObj = writeValue(js)(context, scope)
            ScriptableObject.putProperty(newObj, index, newerObj)

        }
        newObj

      case ESPrimitive(obj) => obj
    }

  /**
    * Runs ECMAScript as JS1.8 (aka ES5)
    *
    * Uses RhinoSandbox to reduce chances injected JS wreaks havoc on the JVM.
    *
    * @see https://en.wikipedia.org/wiki/ECMAScript#Version_correspondence
    */
  def evalRaw(expr: String)(block: (Context, Scriptable) => Unit): AnyRef = {
    // TODO: Parameterize and update callers to pass in source name, max duration, max instructions, etc.?
    // For now, be very liberal with scripts giving 60 seconds of unrestricted CPU usage and unlimited instructions.
    val sourceName = "<ecmascript>"
    val maxDuration: Duration = 60.seconds
    val maxInstructionsOption: Option[Int] = None
    val strict = true
    val languageVersionOption = Option(Context.VERSION_1_8)

    val sandbox = new EnhancedRhinoSandbox(strict, languageVersionOption)
    if (maxDuration.isFinite) {
      sandbox.setMaxDuration(maxDuration.toMillis.toInt)
    }
    maxInstructionsOption foreach sandbox.setInstructionLimit
    sandbox.setUseSafeStandardObjects(true)
    sandbox.setUseSealedScope(true)

    sandbox.eval(sourceName, expr)(block)
  }

  sealed trait ECMAScriptVariable

  case class ESObject(fields: Map[String, ECMAScriptVariable]) extends ECMAScriptVariable
  case class ESArray(array: Array[ECMAScriptVariable]) extends ECMAScriptVariable {
    override def toString: String = s"ESArray(${array.toList})"
  }
  case class ESPrimitive(anyRef: AnyRef) extends ECMAScriptVariable

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
                    encoder: EcmaScriptEncoder,
                    decoder: CwlEcmaScriptDecoder = new CwlEcmaScriptDecoder): ErrorOr[WomValue] = {
    def evaluate = evalRaw(expr) { (context, scope) =>

      val (key, value) = rawValues

      val jsValue = encoder.encode(value)
      val field = writeValue(jsValue)(context, scope)
      ScriptableObject.putProperty(scope, key, field)

      val jsMap = mapValues.mapValues{ _.mapValues(encoder.encode) }

      jsMap foreach {
        case (scopeId, nestedMap) =>

          val newObj = context.newObject(scope)
          nestedMap.toList foreach {
            case (key, value) =>
              val newerObj = writeValue(value)(context, scope)
              ScriptableObject.putProperty(newObj, key, newerObj)
          }
          ScriptableObject.putProperty(scope, scopeId, newObj)
      }

    }
    
    for {
      evaluated <- Try(evaluate).toErrorOr
      decoded <- decoder.decode(evaluated)
    } yield decoded
  }
}
