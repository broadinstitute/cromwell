package cromwell.binding

import cromwell.binding.values._
import cromwell.engine.EngineFunctions
import org.scalatest.mock.MockitoSugar
import org.scalatest.{Matchers, FlatSpec}

import scala.util.{Try, Success, Failure}

class WdlExpressionToFileSpec extends FlatSpec with Matchers with MockitoSugar {
  val expr: String => WdlExpression = WdlExpression.fromString
  def noLookup(String: String): WdlValue = fail("No identifiers should be looked up in this test")

  val engineFunctions = new MockEngineFunctions

  "anonFiles function" should "find no anonymous files in an empty expression" in {
    val e: WdlExpression = expr("")
    e.preevaluateExpressionForFilenames(noLookup, engineFunctions) match {
      case Success(files) => assert(files isEmpty)
      case Failure(error) => fail(error)
    }
  }

  "anonFiles function" should "find no anonymous files in a simple expression" in {
    val e: WdlExpression = expr("1 + 1")
    e.preevaluateExpressionForFilenames(noLookup, engineFunctions) match {
      case Success(files) => assert(files isEmpty)
      case Failure(error) => fail(error)
    }
  }

  "anonFiles function" should "find no anonymous files in stdout() or stderr() expression" in {
    val e: WdlExpression = expr("stdout() + stderr()")
    e.preevaluateExpressionForFilenames(noLookup, engineFunctions) match {
      case Success(files) => assert(files isEmpty)
      case Failure(error) => fail(error)
    }
  }

  "anonFiles function" should "identify read_int's argument as a files" in {
    val e: WdlExpression = expr("""read_int("myfile.txt")""")
    e.preevaluateExpressionForFilenames(noLookup, engineFunctions) match {
      case Success(files) => {
        assert(files.size == 1)
        assert(files.head == WdlFile("myfile.txt"))
      }
      case Failure(error) => fail(error)
    }
  }

  "anonFiles function" should "evaluate expressions within read_int's ( and ) as a single file" in {
    val e: WdlExpression = expr("""read_int("/bin/bash/" + "myfile.txt")""")
    e.preevaluateExpressionForFilenames(noLookup, engineFunctions) match {
      case Success(files) => {
        assert(files.size == 1)
        assert(files.head == WdlFile("/bin/bash/myfile.txt"))
      }
      case Failure(error) => fail(error)
    }
  }

  // I can't see any use for this particular combination but it's possible and should work, so let's test it:
  "anonFiles function" should "evaluate expressions within read_int's which include pre-evaluable function calls" in {
    val e: WdlExpression = expr("""read_int(stdout() + "3")""")
    e.preevaluateExpressionForFilenames(noLookup, engineFunctions) match {
      case Success(files) => {
        assert(files.size == 1)
        assert(files.head == WdlFile("job.stdout.txt/3"))
      }
      case Failure(error) => fail(error)
    }
  }

  "anonFiles function" should "not allow nested read_int functions" in {
    val e: WdlExpression = expr(""" read_int("/etc/" + read_int("somefile") + ".txt")) """)
    e.preevaluateExpressionForFilenames(noLookup, engineFunctions) match {
      case Success(files) => fail("Should not have succeeded with nested read_int functions")
      case Failure(error) => /* expected result; test passes */
    }
  }

  "anonFiles function" should "find file names modified by unary operators" in {
    val e: WdlExpression = expr(""" -read_int("/etc/file1")""")
    e.preevaluateExpressionForFilenames(noLookup, engineFunctions) match {
      case Success(files) => {
        assert(files.size == 1)
        assert(files exists {_.value == "/etc/file1"})
      }
      case Failure(error) => fail(error)
    }
  }

  "anonFiles function" should "find two file names for two consecutive read_X functions" in {
    val e: WdlExpression = expr("""read_int("/etc/file1") + read_string("/bin/file2")""")
    e.preevaluateExpressionForFilenames(noLookup, engineFunctions) match {
      case Success(files) => {
        assert(files.size == 2)
        assert(files exists {_.value == "/etc/file1"})
        assert(files exists {_.value == "/bin/file2"})
      }
      case Failure(error) => fail(error)
    }
  }

  "anonFiles function" should "find all file names in arbitrarily long sequences of read_X functions" in {
    val e: WdlExpression = expr("""read_int("/etc/file1") + read_string("/bin/file2") + read_string("/bin/file3") + read_string("/bin/file4") + read_string("/bin/file5")""")
    e.preevaluateExpressionForFilenames(noLookup, engineFunctions) match {
      case Success(files) => {
        assert(files.size == 5)
        assert(files exists {_.value == "/etc/file1"})
        assert(files exists {_.value == "/bin/file2"})
        assert(files exists {_.value == "/bin/file3"})
        assert(files exists {_.value == "/bin/file4"})
        assert(files exists {_.value == "/bin/file5"})
      }
      case Failure(error) => fail(error)
    }
  }
}

class MockEngineFunctions extends EngineFunctions  {
  // It's important to ensure that these two are never called during pre-evaluation:
  override protected def read_int(params: Seq[Try[WdlValue]]): Try[WdlInteger] = {
    assert(false, "The active 'read_int' function should not be used during pre-evaluation")
    ???
  }
  override protected def read_string(params: Seq[Try[WdlValue]]): Try[WdlString] = {
    assert(false, "The active 'read_string' function should not be used during pre-evaluation")
    ???
  }

  override protected def read_lines(params: Seq[Try[WdlValue]]): Try[WdlArray] = {
    assert(false, "The active 'read_lines' function should not be used during pre-evaluation")
    ???
  }

  // These pre-evaluable functions can be called during the pre-evaluation.
  override protected def stdout(params: Seq[Try[WdlValue]]): Try[WdlFile] = Success(WdlFile("job.stdout.txt"))
  override protected def stderr(params: Seq[Try[WdlValue]]): Try[WdlFile] = Success(WdlFile("job.stderr.txt"))
}
