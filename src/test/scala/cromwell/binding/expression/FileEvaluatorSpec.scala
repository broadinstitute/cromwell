package cromwell.binding.expression

import cromwell.binding.WdlExpression
import cromwell.binding.types._
import cromwell.binding.values._
import org.scalatest.mock.MockitoSugar
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers}

import scala.language.postfixOps

class FileEvaluatorSpec extends FlatSpec with Matchers with MockitoSugar {
  val expr: String => WdlExpression = WdlExpression.fromString
  def noLookup(String: String): WdlValue = fail("No identifiers should be looked up in this test")

  val expressions = Table(
    ("expression", "files", "wdlType"),
    ("1 + 1", Seq.empty[WdlFile], WdlIntegerType),
    ("stdout() + stderr()", Seq.empty[WdlFile], WdlFileType),
    ("""read_int("myfile.txt")""", Seq(WdlSingleFile("myfile.txt")), WdlIntegerType),
    ("""read_int("/bin/bash/" + "myfile.txt")""", Seq(WdlSingleFile("/bin/bash/myfile.txt")), WdlIntegerType),
    ("read_int(stdout())", Seq.empty[WdlFile], WdlIntegerType),
    ("read_int(stdout() + 3)", Seq.empty[WdlFile], WdlIntegerType),
    ("""read_int("/etc/" + read_int("somefile") + ".txt"))""", Seq(WdlSingleFile("somefile")), WdlIntegerType),
    ("""-read_int("/etc/file1")""", Seq(WdlSingleFile("/etc/file1")), WdlIntegerType),
    ("""read_int("/etc/file1") + read_string("/bin/file2")""", Seq(WdlSingleFile("/etc/file1"), WdlSingleFile("/bin/file2")), WdlStringType),
    ("""read_int("/etc/file1") + read_string("/bin/file2") + read_string("/bin/file3") + read_string("/bin/file4") + read_string("/bin/file5")""", Seq(
      WdlSingleFile("/etc/file1"),
      WdlSingleFile("/bin/file2"),
      WdlSingleFile("/bin/file3"),
      WdlSingleFile("/bin/file4"),
      WdlSingleFile("/bin/file5")
    ), WdlStringType),
    (""" "foo" + "bar" """, Seq(WdlSingleFile("foobar")), WdlFileType),
    (""" "foo" + "bar" """, Seq.empty[WdlFile], WdlStringType),
    (""" ["a", "b", "c"] """, Seq(WdlSingleFile("a"), WdlSingleFile("b"), WdlSingleFile("c")), WdlArrayType(WdlFileType)),
    (""" ["a", "b", "c"] """, Seq.empty[WdlFile], WdlArrayType(WdlStringType)),
    (""" {"a": "1", "b": "2", "c":"3"} """, Seq(
      WdlSingleFile("a"),
      WdlSingleFile("1"),
      WdlSingleFile("b"),
      WdlSingleFile("2"),
      WdlSingleFile("c"),
      WdlSingleFile("3")
    ), WdlMapType(WdlFileType, WdlFileType)),
    (""" [read_string("x"), read_string("y")] """, Seq(WdlSingleFile("x"), WdlSingleFile("y")), WdlArrayType(WdlStringType)),
    (""" {read_int("a"): read_string("x"), 4: read_string("y")} """, Seq(WdlSingleFile("a"), WdlSingleFile("x"), WdlSingleFile("y")), WdlArrayType(WdlStringType)),
    (""" glob("out-*.txt") """, Seq(WdlGlobFile("out-*.txt")), WdlFileType)
  )

  forAll (expressions) { (expression, files, wdlType) =>
    it should s"evaluate $expression (coerced to: $wdlType) => $files" in {
      expr(expression).evaluateFiles(noLookup, new NoFunctions, wdlType).get shouldEqual files
    }
  }
}
