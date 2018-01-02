package cwl

import cats.data.NonEmptyList
import wom.types.WomNothingType
import wom.values.{WomArray, WomBoolean, WomFloat, WomInteger, WomMap, WomOptionalValue, WomSingleFile, WomString, WomValue}
import cats.syntax.traverse._
import cats.syntax.apply._
import cats.syntax.validated._
import cats.instances.list._
import common.validation.ErrorOr.ErrorOr
import cats.data.Validated._
import wom.callable.RuntimeEnvironment
import wom.graph.LocalName

import scala.collection.JavaConverters._

object ParameterContext {
  val Empty = ParameterContext()
}

case class ParameterContext(private val inputs: Map[String, AnyRef] = Map.empty,
                       private val self: Array[Map[String, String]] = Array.empty,
                       private val runtime: Map[String, AnyRef] = Map.empty) {

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


  def addInputs(womMap: Map[String, WomValue]): Either[NonEmptyList[String], ParameterContext] = {
    val x: ErrorOr[List[(String, AnyRef)]] = womMap.toList.traverse{
      case (key, womValue) => (key.validNel[String]:ErrorOr[String], toJavascript(womValue)).tupled
    }


    x.map(lst => this.copy(inputs ++ lst.toMap)).toEither
  }

  def addLocalInputs(womMap: Map[LocalName, WomValue]): Either[NonEmptyList[String], ParameterContext] = {
    val stringKeyMap: Map[String, WomValue] = womMap.map {
      case (LocalName(localName), value) => localName -> value
    }
    addInputs(stringKeyMap)
  }

  def ecmaScriptInputs: Map[String, AnyRef] = inputs

  def ecmaScriptCompatinputs: Map[String, AnyRef] = inputs

  def setSelf(newSelf: Array[Map[String, String]]): ParameterContext = this.copy(self = newSelf)

  def ecmaScriptValues:java.util.Map[String, AnyRef] =
    Map(
      "inputs" -> inputs.asJava.asInstanceOf[AnyRef],
      "runtime" -> runtime.asJava.asInstanceOf[AnyRef],
      "self" -> self.asInstanceOf[AnyRef]
    ).asJava

  private def toJavascript(value: WomValue): ErrorOr[AnyRef] = {
    value match {
      case WomOptionalValue(WomNothingType, None) => null
      case WomString(string) => string.validNel
      case WomInteger(int) => int.asInstanceOf[java.lang.Integer].validNel
      case WomFloat(double) => double.asInstanceOf[java.lang.Double].validNel
      case WomBoolean(boolean) => boolean.asInstanceOf[java.lang.Boolean].validNel
      case WomArray(_, array) => array.map(toJavascript).toArray.validNel
      case WomSingleFile(path) => path.validNel
      case WomMap(_, map) =>
        map.toList.traverse({
          case (mapKey, mapValue) => (toJavascript(mapKey), toJavascript(mapValue)).tupled
        }).map(_.toMap.asJava)
      case _ => (s"Value is unsupported in JavaScript: $value").invalidNel
    }
  }
}
