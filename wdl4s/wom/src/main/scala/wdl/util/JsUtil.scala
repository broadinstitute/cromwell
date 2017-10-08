package wdl.util

import javax.script.{ScriptContext, SimpleScriptContext}

import jdk.nashorn.api.scripting.{ClassFilter, NashornScriptEngineFactory, ScriptObjectMirror}
import wdl.types.{WdlArrayType, WdlMapType, WdlNothingType, WdlStringType}
import wdl.values.{WdlArray, WdlBoolean, WdlFloat, WdlInteger, WdlMap, WdlOptionalValue, WdlString, WdlValue}

import scala.collection.JavaConverters._

object JsUtil {

  /**
    * Evaluates a javascript expression.
    *
    * Inputs, and returned output must be one of:
    * - WdlString
    * - WdlBoolean
    * - WdlFloat
    * - WdlInteger
    * - WdlMap
    * - WdlArray
    * - A "WdlNull" equal to WdlOptionalValue(WdlNothingType, None)
    *
    * The WdlMap keys and values, and WdlArray elements must be the one of the above, recursively.
    *
    * WdlSingleFile and WdlGlobFile are not permitted, and must be already converted to one of the above types.
    *
    * @param expr   The javascript expression.
    * @param values A map filled with WDL values.
    * @return The result of the expression.
    */
  def eval(expr: String, values: Map[String, WdlValue] = Map.empty): WdlValue = {
    val engine = ScriptEngineFactory.getScriptEngine(nashornStrictArgs, getNashornClassLoader, noJavaClassFilter)

    val bindings = engine.createBindings()
    val javascriptValues = values.mapValues(toJavascript).asJava
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

  private def toJavascript(value: WdlValue): AnyRef = {
    value match {
      case WdlOptionalValue(WdlNothingType, None) => null
      case WdlString(string) => string
      case WdlInteger(int) => int.asInstanceOf[java.lang.Integer]
      case WdlFloat(double) => double.asInstanceOf[java.lang.Double]
      case WdlBoolean(boolean) => boolean.asInstanceOf[java.lang.Boolean]
      case WdlArray(_, array) => array.map(toJavascript).toArray
      case WdlMap(_, map) =>
        map.map({
          case (mapKey, mapValue) => toJavascript(mapKey) -> toJavascript(mapValue)
        }).asJava
      case _ => throw new IllegalArgumentException(s"Unexpected value: $value")
    }
  }

  private def fromJavascript(value: AnyRef): WdlValue = {
    value match {
      case null => WdlOptionalValue(WdlNothingType, None)
      case string: String => WdlString(string)
      case int: java.lang.Integer => WdlInteger(int)
      case double: java.lang.Double => WdlFloat(double)
      case boolean: java.lang.Boolean => WdlBoolean(boolean)
      case scriptObjectMirror: ScriptObjectMirror if scriptObjectMirror.isArray =>
        val wdlValues = (0 until scriptObjectMirror.size).toArray map { index =>
          fromJavascript(scriptObjectMirror.getSlot(index))
        }
        val wdlArrayType = if (wdlValues.isEmpty) WdlArrayType(WdlNothingType) else WdlArrayType(wdlValues.head.wdlType)
        WdlArray(wdlArrayType, wdlValues)
      case scriptObjectMirror: ScriptObjectMirror if scriptObjectMirror.isFunction =>
        throw new IllegalArgumentException(s"Unexpected function value: $value")
      case scriptObjectMirror: ScriptObjectMirror =>
        val keys = scriptObjectMirror.getOwnKeys(true)
        val wdlMap: Map[WdlString, WdlValue] = keys.map(key =>
          WdlString(key) -> fromJavascript(scriptObjectMirror.get(key))
        ).toMap
        val wdlValues = wdlMap.values
        val wdlValueType = if (wdlValues.isEmpty) WdlArrayType(WdlNothingType) else WdlArrayType(wdlValues.head.wdlType)
        val wdlMapType = WdlMapType(WdlStringType, wdlValueType)
        WdlMap(wdlMapType, keys.map(key =>
          WdlString(key) -> fromJavascript(scriptObjectMirror.get(key))
        ).toMap)
      case _ => throw new IllegalArgumentException(s"Unexpected value: $value")
    }
  }

}
