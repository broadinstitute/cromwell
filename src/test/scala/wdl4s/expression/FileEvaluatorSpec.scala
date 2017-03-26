package wdl4s.expression

import wdl4s.{NoLookup, WdlExpression}
import wdl4s.types._
import wdl4s.values._
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers}

import scala.language.postfixOps

class FileEvaluatorSpec extends FlatSpec with Matchers {
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
    (""" [fileNameAsStringInput, "${fileNameAsStringInput}.bai"] """, Seq(WdlSingleFile("sommat.bam"), WdlSingleFile("sommat.bam.bai")), WdlArrayType(WdlFileType)),
    (""" [ fileNameAsStringInput, mapToFileName["Chris"] ] """, Seq(WdlSingleFile("sommat.bam"), WdlSingleFile("sommatStupid.bam")), WdlArrayType(WdlFileType)),
    (""" {read_int("a"): read_string("x"), 4: read_string("y")} """, Seq(WdlSingleFile("a"), WdlSingleFile("x"), WdlSingleFile("y")), WdlMapType(WdlIntegerType, WdlStringType)),
    (""" glob("out-*.txt") """, Seq(WdlGlobFile("out-*.txt")), WdlFileType),
    (""" glob("out-*.txt")[0] """, Seq(WdlGlobFile("out-*.txt")), WdlFileType),
    (""" read_tsv("my_file") """, Seq(WdlSingleFile("my_file")), WdlFileType),
    (""" if read_int("i.txt") == 10 then "a.txt" else "b.txt" """, Seq(WdlSingleFile("a.txt"), WdlSingleFile("b.txt"), WdlSingleFile("i.txt")), WdlFileType),
    (""" if "a" == "b" then "a.txt" else "b.txt" """, Seq(WdlSingleFile("b.txt")), WdlFileType),
    (""" if b then read_string("t") else "nope" """, Seq(WdlSingleFile("t")), WdlStringType)
  )

  val lookupFunction = Map(
    "b" -> WdlBoolean(true),
    "fileNameAsStringInput" -> WdlString("sommat.bam"),
    "mapToFileName" -> WdlMap(Map(WdlString("Chris") -> WdlString("sommatStupid.bam")))
  )

  forAll (expressions) { (expression, files, wdlType) =>
    it should s"evaluate $expression (coerced to: $wdlType) => $files" in {
      WdlExpression.fromString(expression).evaluateFiles(lookupFunction, NoFunctions, wdlType).get.toSet should be(files.toSet)
    }
  }
}
