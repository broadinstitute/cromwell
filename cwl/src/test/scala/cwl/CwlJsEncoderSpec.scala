package cwl

import common.validation.Validation._
import cwl.internal.JsUtil
import cwl.internal.JsUtil.{JsArray, JsObject, JsPrimitive}
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
      "class" -> JsPrimitive("File"),
      "location" -> JsPrimitive("path/to/file.txt"),
      "path" -> JsPrimitive("path/to/file.txt"),
      "dirname" -> JsPrimitive("path/to"),
      "basename" -> JsPrimitive("file.txt"),
      "nameroot" -> JsPrimitive("file"),
      "nameext" -> JsPrimitive(".txt")
    )
    val result: JsUtil.Js = encoder.encode(file)
    val resultMap = result.asInstanceOf[JsObject].fields
    resultMap.filterKeys(_ != "secondaryFiles") should contain theSameElementsAs expected
    resultMap("secondaryFiles") should be(a[JsArray])
    resultMap("secondaryFiles").asInstanceOf[JsArray].array should be(empty)
  }

}
