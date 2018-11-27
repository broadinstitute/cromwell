package cwl

import common.validation.Validation._
import cwl.internal.{EcmaScriptEncoder, EcmaScriptUtil}
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wom.values.{WomFloat, WomMaybePopulatedFile, WomString, WomValue}

class CwlEcmaScriptDecoderSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
  behavior of "CwlJsDecoder"

  val decodeTests= Table[String, String, Map[String, Map[String, WomValue]], WomValue](
    (
      "description",
      "expr",
      "values",
      "expected"
    ),
    (
      "a passed through file",
      "inputMap",
      Map("inputMap" -> Map(
        "class" -> WomString("File"),
        "location" -> WomString("path/to/file.txt"),
        "checksum" -> WomString("hash_here"),
        "size" -> WomFloat(Double.box(123.4)),
        "format" -> WomString("file_format"),
        "contents" -> WomString("file_contents")
      )),
      WomMaybePopulatedFile(
        valueOption = Option("path/to/file.txt"),
        checksumOption = Option("hash_here"),
        sizeOption = Option(123L),
        formatOption = Option("file_format"),
        contentsOption = Option("file_contents"),
        secondaryFiles = Vector()
      )
    ),
    (
      "a script generated file",
      """|var result = {
         |"class": "File",
         |"path": "other/dir/img.jpg",
         |"checksum": "check_sum",
         |"size": 567.8,
         |"format": "more_jpeg",
         |"contents": "lol_cat"
         |};
         |result;
         |""".stripMargin,
      Map(),
      WomMaybePopulatedFile(
        valueOption = Option("other/dir/img.jpg"),
        checksumOption = Option("check_sum"),
        sizeOption = Option(567L),
        formatOption = Option("more_jpeg"),
        contentsOption = Option("lol_cat"),
        secondaryFiles = Vector()
      )
    )
  )

  forAll(decodeTests) { (description, expr, values: Map[String, Map[String, WomValue]], expected) =>
    it should s"decode $description" in {
      val result = EcmaScriptUtil.evalStructish(
        expr,
        "fake" -> WomString("unused"),
        values,
        new EcmaScriptEncoder()
      ).toTry.get
      result should be(expected)
    }
  }
}
