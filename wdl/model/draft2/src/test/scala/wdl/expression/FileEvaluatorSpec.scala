package wdl.expression

import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.{FlatSpec, Matchers}
import wdl.WdlExpression
import wom.types._
import wom.values._

class FileEvaluatorSpec extends FlatSpec with Matchers {
  val expressions = Table(
    ("expression", "files", "womType"),
    ("1 + 1", Seq.empty[WomFile], WomIntegerType),
    ("stdout() + stderr()", Seq.empty[WomFile], WomSingleFileType),
    ("""read_int("myfile.txt")""", Seq(WomSingleFile("myfile.txt")), WomIntegerType),
    ("""read_int("/bin/bash/" + "myfile.txt")""", Seq(WomSingleFile("/bin/bash/myfile.txt")), WomIntegerType),
    ("read_int(stdout())", Seq.empty[WomSingleFile], WomIntegerType),
    ("read_int(stdout() + 3)", Seq.empty[WomSingleFile], WomIntegerType),
    ("""read_int("/etc/" + read_int("somefile") + ".txt"))""", Seq(WomSingleFile("somefile")), WomIntegerType),
    ("""-read_int("/etc/file1")""", Seq(WomSingleFile("/etc/file1")), WomIntegerType),
    ("""read_int("/etc/file1") + read_string("/bin/file2")""", Seq(WomSingleFile("/etc/file1"), WomSingleFile("/bin/file2")), WomStringType),
    ("""read_int("/etc/file1") + read_string("/bin/file2") + read_string("/bin/file3") + read_string("/bin/file4") + read_string("/bin/file5")""", Seq(
      WomSingleFile("/etc/file1"),
      WomSingleFile("/bin/file2"),
      WomSingleFile("/bin/file3"),
      WomSingleFile("/bin/file4"),
      WomSingleFile("/bin/file5")
    ), WomStringType),
    (""" "foo" + "bar" """, Seq(WomSingleFile("foobar")), WomSingleFileType),
    (""" "foo" + "bar" """, Seq.empty[WomSingleFile], WomStringType),
    (""" ["a", "b", "c"] """,
      Seq(WomSingleFile("a"), WomSingleFile("b"), WomSingleFile("c")), WomArrayType(WomSingleFileType)),
    (""" ["a", "b", "c"] """, Seq.empty[WomSingleFile], WomArrayType(WomStringType)),
    (""" {"a": "1", "b": "2", "c":"3"} """, Seq(
      WomSingleFile("a"),
      WomSingleFile("1"),
      WomSingleFile("b"),
      WomSingleFile("2"),
      WomSingleFile("c"),
      WomSingleFile("3")
    ), WomMapType(WomSingleFileType, WomSingleFileType)),
    (""" [read_string("x"), read_string("y")] """, Seq(WomSingleFile("x"), WomSingleFile("y")), WomArrayType(WomStringType)),
    (""" [fileNameAsStringInput, "${fileNameAsStringInput}.bai"] """,
      Seq(WomSingleFile("sommat.bam"), WomSingleFile("sommat.bam.bai")), WomArrayType(WomSingleFileType)),
    (""" [ fileNameAsStringInput, mapToFileName["Chris"] ] """,
      Seq(WomSingleFile("sommat.bam"), WomSingleFile("sommatStupid.bam")), WomArrayType(WomSingleFileType)),
    (""" {read_int("a"): read_string("x"), 4: read_string("y")} """,
      Seq(WomSingleFile("a"), WomSingleFile("x"), WomSingleFile("y")), WomMapType(WomIntegerType, WomStringType)),
    (""" glob("out-*.txt") """, Seq(WomGlobFile("out-*.txt")), WomSingleFileType),
    (""" glob("out-*.txt")[0] """, Seq(WomGlobFile("out-*.txt")), WomSingleFileType),
    (""" read_tsv("my_file") """, Seq(WomSingleFile("my_file")), WomSingleFileType),
    (""" if read_int("i.txt") == 10 then "a.txt" else "b.txt" """,
      Seq(WomSingleFile("a.txt"), WomSingleFile("b.txt"), WomSingleFile("i.txt")), WomSingleFileType),
    (""" if "a" == "b" then "a.txt" else "b.txt" """, Seq(WomSingleFile("b.txt")), WomSingleFileType),
    (""" if b then read_string("t") else "nope" """, Seq(WomSingleFile("t")), WomStringType),
    (""" read_string(basename(fileInput, ".txt") + ".bam") """, Seq(WomSingleFile("input.bam")), WomStringType),
    (""" size("foo.txt") """, Seq(WomSingleFile("foo.txt")), WomFloatType),
    (""" round(size("foo.txt")) """, Seq(WomSingleFile("foo.txt")), WomIntegerType),
    (""" size("foo.txt", "GB") """, Seq(WomSingleFile("foo.txt")), WomIntegerType),
    (""" round(size("foo.txt", "GB")) """, Seq(WomSingleFile("foo.txt")), WomIntegerType)

  )

  val lookupFunction = Map(
    "b" -> WomBoolean(true),
    "fileInput" -> WomSingleFile("gs://bucket/path/to/input.txt"),
    "fileNameAsStringInput" -> WomString("sommat.bam"),
    "mapToFileName" -> WomMap(Map(WomString("Chris") -> WomString("sommatStupid.bam")))
  )

  forAll (expressions) { (expression, files, anticipatedType) =>
    it should s"evaluate $expression (coerced to: $anticipatedType) => $files" in {
      WdlExpression.fromString(expression).evaluateFiles(lookupFunction, PureStandardLibraryFunctions, anticipatedType).get.toSet should be(files.toSet)
    }
  }
}
