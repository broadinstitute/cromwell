package cwl

import wom.callable.RuntimeEnvironment
import wom.graph.LocalName
import wom.types.WomNothingType
import wom.values.{WomArray, WomBoolean, WomFloat, WomInteger, WomMap, WomOptionalValue, WomSingleFile, WomString, WomValue}

import scala.collection.JavaConverters._

object ParameterContext {

  val Empty = ParameterContext()

  def apply(inputs: Map[String, WomValue]): ParameterContext = ParameterContext().addInputs(inputs)
}

/**
  * This class holds values that ultimately will be bound to variables in JDK's "Nashorn" Javascript execution engine.
  *
  * The untyped format accomplishes two goals:
  *   * Allows multiple types to be stored in one map
  *   * Holds values as Nashorn expects them (a map)
  *
  * See {{cwl.JsUtil}} for nashorn-specific implementation details.
  *
  */
case class ParameterContext(private val inputs: Map[String, AnyRef] = Map.empty,
                            private val self: Array[Map[String, String]] = Array.empty,
                            private val runtime: Map[String, AnyRef] = Map.empty) {

  /**
    * see <a href="http://www.commonwl.org/v1.0/CommandLineTool.html#Runtime_environment">CWL Spec</a>.
    */
  def setRuntime(runtimeEnvironment: RuntimeEnvironment): ParameterContext = {
    import runtimeEnvironment._

    this.copy(runtime =  Map(
      "outdir" -> outputPath,
      "tmpdir" -> tempPath,
      "cores" -> cores.toString,
      "ram" -> ram.toString,
      "outdirSize" -> outputPathSize.toString,
      "tmpdirSize" -> tempPathSize.toString
    ))
  }

  /**
    * The public method of adding values hides our untyped internal representation and verifies that a Nashorn-compatible
    * value exists for the womvalue being passed in.
    *
    * @param womMap
    * @return a new Parameter Context with the inputs values added
    */
  def addInputs(womMap: Map[String, WomValue]): ParameterContext = {
    val jsValues = womMap.mapValues(toJavascript)
    this.copy(inputs ++ jsValues.toMap)
  }

  def addLocalInputs(womMap: Map[LocalName, WomValue]): ParameterContext = {
    val stringKeyMap: Map[String, WomValue] = womMap.map {
      case (LocalName(localName), value) => localName -> value
    }
    addInputs(stringKeyMap)
  }

  /**
    * Expose the inputs for logging convenience.
    */
  def ecmaScriptInputs: Map[String, AnyRef] = inputs

  //TODO CWL: This will change as Self is a shapeshifter.  See http://www.commonwl.org/v1.0/CommandLineTool.html#Parameter_references
  def setSelf(newSelf: Array[Map[String, String]]): ParameterContext = this.copy(self = newSelf)

  /**
    * Outputs the values in a java map as Nashorn expects it with keys as CWL expects them.
    *
    * For the static key values see http://www.commonwl.org/v1.0/CommandLineTool.html#Parameter_references
    */
  def ecmaScriptValues:java.util.Map[String, AnyRef] =
    Map(
      "inputs" -> inputs.asJava.asInstanceOf[AnyRef],
      "runtime" -> runtime.asJava.asInstanceOf[AnyRef],
      "self" -> self.map(_.asJava).asInstanceOf[AnyRef]
    ).asJava

  //Nashorn expects values in Java native format
  private def toJavascript(value: WomValue): AnyRef = {
    value match {
      case WomOptionalValue(WomNothingType, None) => null
      case WomString(string) => string
      case WomInteger(int) => Int.box(int)
      case WomFloat(double) => Double.box(double)
      case WomBoolean(boolean) => Boolean.box(boolean)
      case WomArray(_, array) => array.map(toJavascript).toArray
      case WomSingleFile(path) => path
      case WomMap(_, map) =>
        map.map{
          case (mapKey, mapValue) => toJavascript(mapKey) -> toJavascript(mapValue)
        }.toMap.asJava
      case _ => throw new RuntimeException(s"Value is unsupported in JavaScript: $value")
    }
  }
}
