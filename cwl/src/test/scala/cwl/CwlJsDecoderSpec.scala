package cwl

import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import common.validation.Validation._
import wom.util.JsUtil
import wom.values.{WomMaybePopulatedFile, WomValue}

import scala.collection.JavaConverters._

class CwlJsDecoderSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
  behavior of "CwlJsDecoder"

  val decodeTests = Table[String, String, Map[String, AnyRef], WomValue](
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
        "class" -> "File",
        "location" -> "path/to/file.txt",
        "checksum" -> "hash_here",
        "size" -> Double.box(123.4),
        "format" -> "file_format",
        "contents" -> "file_contents"
      ).asJava),
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

  forAll(decodeTests) { (description, expr, values, expected) =>
    it should s"decode $description" in {
      val decoder = new CwlJsDecoder
      val rawResult = JsUtil.evalRaw(expr, values.asJava).toTry.get
      val result = decoder.decode(rawResult).toTry.get
      result should be(expected)
    }
  }
}
