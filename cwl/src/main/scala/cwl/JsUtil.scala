package cwl

import javax.script.{ScriptContext, SimpleScriptContext}

import jdk.nashorn.api.scripting.{ClassFilter, NashornScriptEngineFactory, ScriptObjectMirror}
import wom.types._
import wom.values._

import scala.collection.JavaConverters._

object JsUtil {

  /**
    * Evaluates a javascript expression.
    *
    * Inputs, and returned output must be one of:
    * - WomString
    * - WomBoolean
    * - WomFloat
    * - WomInteger
    * - WomMap
    * - WomArray
    * - A "WomNull" equal to WomOptionalValue(WomNothingType, None)
    *
    * The WomMap keys and values, and WomArray elements must be the one of the above, recursively.
    *
    * WomSingleFile and WomGlobFile are not permitted, and must be already converted to one of the above types.
    *
    * @param expr   The javascript expression.
    * @param values A map filled with WOM values.
    * @return The result of the expression.
    */
  def eval(expr: String, javascriptValues: java.util.Map[String, AnyRef]): WomValue = {
    println(s"evaluating $expr against $javascriptValues")
    val engine = ScriptEngineFactory.getScriptEngine(nashornStrictArgs, getNashornClassLoader, noJavaClassFilter)

    val bindings = engine.createBindings()
    bindings.asInstanceOf[java.util.Map[String, Any]].putAll(javascriptValues)

    val context = new SimpleScriptContext
    context.setBindings(bindings, ScriptContext.ENGINE_SCOPE)
    val result = engine.eval(expr, context)
    fromJavascript(result)
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
  private val nashornStrictArgs = Array("-doe", "-strict", "--no-java", "--no-syntax-extensions", "--language=es5")

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

  private def toJavascript(value: WomValue): AnyRef = {
    value match {
      case WomOptionalValue(WomNothingType, None) => null
      case WomString(string) => string
      case WomInteger(int) => int.asInstanceOf[java.lang.Integer]
      case WomFloat(double) => double.asInstanceOf[java.lang.Double]
      case WomBoolean(boolean) => boolean.asInstanceOf[java.lang.Boolean]
      case WomArray(_, array) => array.map(toJavascript).toArray
      case WomSingleFile(path) => path
      case WomMap(_, map) =>
        map.map({
          case (mapKey, mapValue) => toJavascript(mapKey) -> toJavascript(mapValue)
        }).asJava
      case _ => throw new IllegalArgumentException(s"Unexpected value: $value")
    }
  }

  private def fromJavascript(value: AnyRef): WomValue = {
    println(s" from javascript $value")
    def isWhole(d: Double) = (d == Math.floor(d)) && !java.lang.Double.isInfinite(d)
    value match {
      case null => WomOptionalValue(WomNothingType, None)
      case string: String => WomString(string)
      case int: java.lang.Integer => WomInteger(int)
      case int: java.lang.Double if isWhole(int) => WomInteger(int.intValue()) // Because numbers in nashorn come back as 'Double's
      case double: java.lang.Double => WomFloat(double)
      case boolean: java.lang.Boolean => WomBoolean(boolean)
      case scriptObjectMirror: ScriptObjectMirror if scriptObjectMirror.isArray =>
        val womValues = (0 until scriptObjectMirror.size).toArray map { index =>
          fromJavascript(scriptObjectMirror.getSlot(index))
        }
        val womArrayType = if (womValues.isEmpty) WomArrayType(WomNothingType) else WomArrayType(womValues.head.womType)
        WomArray(womArrayType, womValues)
      case scriptObjectMirror: ScriptObjectMirror if scriptObjectMirror.isFunction =>
        throw new IllegalArgumentException(s"Unexpected function value: $value")
      case scriptObjectMirror: ScriptObjectMirror =>
        val keys = scriptObjectMirror.getOwnKeys(true)
        val womMap: Map[WomString, WomValue] = keys.map(key =>
          WomString(key) -> fromJavascript(scriptObjectMirror.get(key))
        ).toMap
        val womValues = womMap.values
        val womValueType = if (womValues.isEmpty) WomArrayType(WomNothingType) else WomArrayType(womValues.head.womType)
        val womMapType = WomMapType(WomStringType, womValueType)
        WomMap(womMapType, keys.map(key =>
          WomString(key) -> fromJavascript(scriptObjectMirror.get(key))
        ).toMap)
      case _ => throw new IllegalArgumentException(s"Unexpected value: $value")
    }
  }
}
