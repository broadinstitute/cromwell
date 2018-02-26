package cwl

import common.validation.Validation._
import cwl.internal.JsUtil
import cwl.internal.JsUtil.{ESArray, ESObject, ESPrimitive}
import org.mozilla.javascript.NativeObject
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wom.values.WomMaybePopulatedFile

import scala.collection.JavaConverters._

class CwlJsEncoderSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "CwlJsEncoder"

  it should "encode" in {
    val encoder = new CwlJsEncoder
    val file = WomMaybePopulatedFile("path/to/file.txt")
    val expected = Map(
      "class" -> ESPrimitive("File"),
      "location" -> ESPrimitive("path/to/file.txt"),
      "path" -> ESPrimitive("path/to/file.txt"),
      "dirname" -> ESPrimitive("path/to"),
      "basename" -> ESPrimitive("file.txt"),
      "nameroot" -> ESPrimitive("file"),
      "nameext" -> ESPrimitive(".txt")
    )
    val result: JsUtil.ECMAScriptVariable = encoder.encode(file)
    val resultMap = result.asInstanceOf[ESObject].fields
    resultMap.filterKeys(_ != "secondaryFiles") should contain theSameElementsAs expected
    resultMap("secondaryFiles") should be(a[ESArray])
    resultMap("secondaryFiles").asInstanceOf[ESArray].array should be(empty)
  }

}
