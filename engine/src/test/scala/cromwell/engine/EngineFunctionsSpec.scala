package cromwell.engine

import java.nio.file.Path

import cromwell.backend.wdl.{ReadLikeFunctions, WriteFunctions}
import cromwell.core.path.{DefaultPathBuilder, PathBuilder}
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.Tables.Table
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.expression.{NoFunctions, PureStandardLibraryFunctionsLike, WdlStandardLibraryFunctions}
import wdl4s.values.{WdlFile, WdlInteger, WdlString, WdlValue}

import scala.util.{Failure, Success, Try}

class EngineFunctionsSpec extends FlatSpec with Matchers {

  trait WdlStandardLibraryImpl extends WdlStandardLibraryFunctions with ReadLikeFunctions with WriteFunctions with PureStandardLibraryFunctionsLike {
    private def fail(name: String) = Failure(new NotImplementedError(s"$name() not implemented yet"))

    override def writeTempFile(path: String, prefix: String, suffix: String, content: String): String = super[WriteFunctions].writeTempFile(path, prefix, suffix, content)
    override def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stdout")
    override def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = fail("stderr")
  }

  def expectFailure(value: Try[WdlValue]) = value match {
    case Success(s) => fail(s"$s: Expected this function invocation to fail")
    case Failure(ex) => // expected
  }
  "EngineFunctions" should "all initially be undefined" in {
    val stdFunctions = Seq(
      "stdout", "stderr", "read_lines", "read_tsv", "read_map", "read_object", "read_objects",
      "read_json", "read_int", "read_string", "read_float", "read_boolean", "write_lines",
      "write_tsv", "write_map", "write_object", "write_objects", "write_json", "size", "sub"
    )
    stdFunctions.foreach {func =>
      expectFailure(NoFunctions.getFunction(func)(Seq.empty[Try[WdlValue]]))
    }
  }

  "sub" should "replace a string according to a pattern" in {
    class TestEngineFn extends WdlStandardLibraryImpl {
      override def glob(path: String, pattern: String): Seq[String] = ???
      override def pathBuilders: List[PathBuilder] = List(DefaultPathBuilder)
      override def writeDirectory: Path = ???
    }

    val engineFn = new TestEngineFn

    val table = Table(
      ("str", "pattern", "replace", "result"),
      ("somefilename.bam", ".bam$", ".txt", "somefilename.txt"),
      ("somefilename.bam", "\\.bam$", "", "somefilename"),
      ("somefilename.bam.bam", ".bam$", ".txt", "somefilename.bam.txt"),
      ("somefilenamewith.baminside.ext", ".bam$", ".txt", "somefilenamewith.baminside.ext"),
      ("gs://some/gcs/path/to/my_bame_file.bam", "gs://.*/", "", "my_bame_file.bam"),
      ("somefilename", "^some", "other", "otherfilename")
    )

    forAll(table) { (str, pattern, replace, result) =>
      val stringSubstitution: Try[WdlString] = engineFn.sub(List(Success(WdlString(str)), Success(WdlString(pattern)), Success(WdlString(replace))))
      stringSubstitution.isSuccess shouldBe true
      stringSubstitution.get.valueString shouldBe result

      val fileSubstitution: Try[WdlString] = engineFn.sub(List(Success(WdlFile(str)), Success(WdlString(pattern)), Success(WdlString(replace))))
      fileSubstitution.isSuccess shouldBe true
      fileSubstitution.get.valueString shouldBe result
    }

    val nonCompliantValues = List(
      Seq(Success(WdlFile("input")), Success(WdlInteger(1)), Success(WdlString("replaces"))),
      Seq(Failure(new Exception), Success(WdlString("pattern")), Success(WdlString("replaces")))
    )

    nonCompliantValues foreach { nonCompliantValue =>
      val sub: Try[WdlString] = engineFn.sub(nonCompliantValue)
      sub.isFailure shouldBe true
      val failure = sub.failed.get
      failure.isInstanceOf[IllegalArgumentException] shouldBe true
      failure.getMessage should include("Invalid parameters for engine function sub")
    }

    val failedSub = engineFn.sub(nonCompliantValues.head :+ Success(WdlString("extra value")))
    failedSub.isFailure shouldBe true
    val failure = failedSub.failed.get
    failure.isInstanceOf[IllegalArgumentException] shouldBe true
    failure.getMessage shouldBe "Invalid number of parameters for engine function sub: 4. sub takes exactly 3 parameters."
  }
}
