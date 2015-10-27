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
    ("""read_int("myfile.txt")""", Seq(WdlFile("myfile.txt")), WdlIntegerType),
    ("""read_int("/bin/bash/" + "myfile.txt")""", Seq(WdlFile("/bin/bash/myfile.txt")), WdlIntegerType),
    ("read_int(stdout())", Seq.empty[WdlFile], WdlIntegerType),
    ("read_int(stdout() + 3)", Seq.empty[WdlFile], WdlIntegerType),
    ("""read_int("/etc/" + read_int("somefile") + ".txt"))""", Seq(WdlFile("somefile")), WdlIntegerType),
    ("""-read_int("/etc/file1")""", Seq(WdlFile("/etc/file1")), WdlIntegerType),
    ("""read_int("/etc/file1") + read_string("/bin/file2")""", Seq(WdlFile("/etc/file1"), WdlFile("/bin/file2")), WdlStringType),
    ("""read_int("/etc/file1") + read_string("/bin/file2") + read_string("/bin/file3") + read_string("/bin/file4") + read_string("/bin/file5")""", Seq(
      WdlFile("/etc/file1"),
      WdlFile("/bin/file2"),
      WdlFile("/bin/file3"),
      WdlFile("/bin/file4"),
      WdlFile("/bin/file5")
    ), WdlStringType),
    (""" "foo" + "bar" """, Seq(WdlFile("foobar")), WdlFileType),
    (""" "foo" + "bar" """, Seq.empty[WdlFile], WdlStringType),
    (""" ["a", "b", "c"] """, Seq(WdlFile("a"), WdlFile("b"), WdlFile("c")), WdlArrayType(WdlFileType)),
    (""" ["a", "b", "c"] """, Seq.empty[WdlFile], WdlArrayType(WdlStringType)),
    (""" {"a": "1", "b": "2", "c":"3"} """, Seq(
      WdlFile("a"),
      WdlFile("1"),
      WdlFile("b"),
      WdlFile("2"),
      WdlFile("c"),
      WdlFile("3")
    ), WdlMapType(WdlFileType, WdlFileType)),
    (""" [read_string("x"), read_string("y")] """, Seq(WdlFile("x"), WdlFile("y")), WdlArrayType(WdlStringType)),
    (""" {read_int("a"): read_string("x"), 4: read_string("y")} """, Seq(WdlFile("a"), WdlFile("x"), WdlFile("y")), WdlArrayType(WdlStringType))
  )

  forAll (expressions) { (expression, files, wdlType) =>
    it should s"evaluate $expression (coerced to: $wdlType) => $files" in {
      expr(expression).evaluateFiles(noLookup, new NoFunctions, wdlType).get shouldEqual files
    }
  }
}
