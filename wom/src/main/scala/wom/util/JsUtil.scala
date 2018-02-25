package wom.util

import javax.script.{ScriptContext, SimpleScriptContext}

import common.validation.ErrorOr._
import common.validation.Validation._
import jdk.nashorn.api.scripting.{ClassFilter, NashornScriptEngineFactory}
import wom.values._

import scala.collection.JavaConverters._

object JsUtil {
  /**
    * Evaluates a javascript expression using WOM values.
    *
    * @param expr    The javascript expression.
    * @param values  A map filled with WOM values.
    * @param encoder Encodes wom values to javascript.
    * @param decoder Decodes wom values from javascript.
    * @return The result of the expression.
    */
  def eval(expr: String,
           values: Map[String, WomValue] = Map.empty,
           encoder: JsEncoder = new JsEncoder,
           decoder: JsDecoder = new JsDecoder): ErrorOr[WomValue] = {
    for {
      javaScriptValues <- values.traverseValues(encoder.encode)
      result <- evalRaw(expr, javaScriptValues.asJava)
      decoded <- decoder.decode(result)
    } yield decoded
  }

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
                    rawValues: Map[String, WomValue] = Map.empty,
                    mapValues: Map[String, Map[String, WomValue]] = Map.empty,
                    encoder: JsEncoder = new JsEncoder,
                    decoder: JsDecoder = new JsDecoder): ErrorOr[WomValue] = {
    for {
      rawJavaScriptValues <- rawValues.traverseValues(encoder.encode)
      mapJavaScriptValues <- mapValues
        .mapValues(_.traverseValues(encoder.encode).map(JsMap))
        .traverseValues(identity)
      javaScriptValues = rawJavaScriptValues ++ mapJavaScriptValues
      result <- evalRaw(expr, javaScriptValues.asJava)
      decoded <- decoder.decode(result)
    } yield decoded
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
