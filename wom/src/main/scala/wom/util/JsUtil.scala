package wom.util

import javax.script.{ScriptContext, SimpleScriptContext}

import cats.data.{Kleisli, NonEmptyList, Validated}
import common.validation.ErrorOr._
import common.validation.Validation._
import jdk.nashorn.api.scripting.{ClassFilter, NashornScriptEngineFactory}
import org.mozilla.javascript.{Context, Scriptable, ScriptableObject}
import wom.values._
import cats.syntax.traverse._
import cats.instances.list._
import cats.data.Kleisli._
import cats.data.Validated._
import cats.data.Validated._

import scala.collection.JavaConverters._

object JsUtil {
  val encoder = new JsEncoder
  def writeValue(value: Js)(scriptable: Scriptable, context: Context, scope: ScriptableObject): AnyRef =
    value match {
      case JsObject(fields) => {
        val newObj = context.newObject(scope)

        fields.toList.foreach{
          case (name, value) =>
            val newerObj = writeValue(value)(newObj, context, scope)
            ScriptableObject.putProperty(newObj, name, newerObj)
        }
        newObj
      }

      case JsArray(name, array) =>
        val newObj = context.newArray(scope, array.length)

        array.toList.zipWithIndex.foreach {
          case (js, index) =>

            val newerObj = writeValue(js)
            ScriptableObject.putProperty(newObj, index, newerObj)

        }
        ScriptableObject.putProperty(scope, name, newObj)
        newObj

      case JsField(obj) =>
        obj
    }

  def writeValue(name: String, value: WomValue)(scriptable: Scriptable, context: Context, scope: ScriptableObject): Validated[NonEmptyList[String], Unit] =
    encoder.encode(value).map{writeValue(_)(scriptable, context, scope) }

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
  case class JsArray(name: String, array: Array[Js]) extends Js
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
    evalRaw(expr) { (context, so) =>

      val obj = context.newObject(so)


    for {
      rawJavaScriptValues <- rawValues.toList.traverse[ErrorOr, Unit]{
        case (key, value) => writeValue(key, value)(obj)
      }
      mapJavaScriptValues <- mapValues.mapValues(_.traverseValues(encoder.encode).map(JsMap)).traverseValues(identity)
      javaScriptValues = rawJavaScriptValues ++ mapJavaScriptValues
      result <- evalRaw(expr, javaScriptValues.asJava)
      decoded <- decoder.decode(result)
    } yield decoded
    }
  }

  /**
    * Evaluates a javascript expression using raw javascript compatible values.
    *
    * @param expr   The javascript expression.
    * @param values A map filled with javascript compatible values.
    * @return The result of the expression.
    */
  def evalRaw(expr: String, values: java.util.Map[String, AnyRef]): ErrorOr[AnyRef] = {
    validate {
      val engine = ScriptEngineFactory.getScriptEngine(nashornStrictArgs, getNashornClassLoader, noJavaClassFilter)
      val bindings = engine.createBindings()
      bindings.asInstanceOf[java.util.Map[String, Any]].putAll(values)

      val context = new SimpleScriptContext
      context.setBindings(bindings, ScriptContext.ENGINE_SCOPE)
      engine.eval(expr, context)
    }
  }

  private val ScriptEngineFactory = new NashornScriptEngineFactory

  /**
    * Duplicate of the private jdk.nashorn.api.scripting.NashornScriptEngineFactory#getAppClassLoader()
    */
  private def getNashornClassLoader = {
    Option(Thread.currentThread.getContextClassLoader).getOrElse(classOf[NashornScriptEngineFactory].getClassLoader)
  }

  /**
    * Add stricter nashorn arguments to the default `-doe`.
    *
    * https://docs.oracle.com/javase/8/docs/technotes/tools/windows/jjs.html
    * https://github.com/JetBrains/jdk8u_nashorn/blob/jdk8u76-b03/src/jdk/nashorn/internal/runtime/resources/Options.properties
    * [[jdk.nashorn.api.scripting.NashornScriptEngineFactory#DEFAULT_OPTIONS]]
    */
  private val nashornStrictArgs = Array("-doe", "--no-java", "-strict", "--no-syntax-extensions", "--language=es5")

  /**
    * Don't allow any java classes.
    */
  private val noJavaClassFilter = {
    // TODO: when 2.12+...
    //private val noJavaClassFilter: ClassFilter = _ => false
    new ClassFilter {
      override def exposeToScripts(unused: String): Boolean = false
    }
  }
}
