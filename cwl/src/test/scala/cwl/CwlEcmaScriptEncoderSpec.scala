package cwl

import cwl.internal.EcmaScriptUtil.{ESArray, ESObject, ESPrimitive}
import cwl.internal.{EcmaScriptEncoder, EcmaScriptUtil}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wom.values.WomMaybePopulatedFile

class CwlEcmaScriptEncoderSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "EcmaScriptEncoder"

  it should "encode" in {
    val encoder = new EcmaScriptEncoder()
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
    val result: EcmaScriptUtil.ECMAScriptVariable = encoder.encode(file)
    val resultMap = result.asInstanceOf[ESObject].fields
    resultMap.filterKeys(_ != "secondaryFiles").toList should contain theSameElementsAs expected
    resultMap("secondaryFiles") should be(a[ESArray])
    resultMap("secondaryFiles").asInstanceOf[ESArray].array should be(empty)
  }

}
