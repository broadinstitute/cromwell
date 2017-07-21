package cromwell.backend.wdl

import java.nio.file.{Paths, Path}
import java.nio.file.StandardOpenOption._
import scala.util.{Failure, Success, Try}

import cromwell.backend.standard.{DefaultStandardExpressionFunctionsParams, StandardExpressionFunctions}
import cromwell.core.CallContext
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.wdl.values._
import com.google.common.io.Files
import fs2.{Task, Stream}

class FileSizeSpec extends FlatSpec with Matchers {
  val _readLinesLimit = 4
  val _readBoolLimit = 5
  val _readIntLimit = 6
  val _readFloatLimit = 7
  val _readStringLimit = 8
  val _readJsonLimit = 9
  val _readTsvLimit = 10
  val _readMapLimit = 11
  val _readObjectLimit = 12

  val rlf = {
    val path = DefaultPathBuilder.build("/tmp").get

    val dp = DefaultStandardExpressionFunctionsParams(List(cromwell.core.path.DefaultPathBuilder), CallContext(path, "stdout", "stderr"))

    new StandardExpressionFunctions(dp) {
      override val fileSizeLimitationConfig =
        new FileSizeLimitationConfig {
          val readLinesLimit = _readLinesLimit
          val readIntLimit = _readIntLimit
          val readFloatLimit = _readFloatLimit
          val readStringLimit = _readStringLimit
          val readJsonLimit = _readJsonLimit
          val readTsvLimit = _readTsvLimit
          val readMapLimit = _readMapLimit
          val readObjectLimit = _readObjectLimit
          val readBoolLimit = _readBoolLimit
        }
    }
  }

  val tempDir = Files.createTempDir
  tempDir.deleteOnExit

  def testOverUnder(command: String, n: Int, f: ReadLikeFunctions => (Seq[Try[WdlValue]] => Try[WdlValue])) = {

    def testInner(n: Int, test: Try[WdlValue] => Unit) = {

      def createTempFileOfSize(size: Int): Path = {

        val fn = tempDir.toString + "/" + scala.util.Random.alphanumeric.take(5).mkString
        val jPath = Paths.get(fn)
        jPath.toFile.deleteOnExit
        val start = Stream[Task, Byte](1).repeat.take(size.toLong)
        val end = fs2.io.file.writeAll[Task](jPath, Seq(CREATE_NEW, WRITE))
        (start to end).run.unsafeRunSync
        //jPath is now a file of n bytes, we can return it
        jPath
      }

      val file = createTempFileOfSize(n)
      val params = Seq(Try(WdlString(file.toString)))

      f(rlf)(params) match {
        case t => test(t)
      }
    }

    def testOver() = {
      testInner(n + 1, {
        case Failure(_: FileSizeTooBig) => //success
        case t => throw new RuntimeException(s"should not have eaten this file that is too big! msg: $t")
      })
    }

    def testUnder() = {
      testInner(n - 1, {
        case Success(_) =>
        case Failure(_: NumberFormatException) => //we're not testing parsing
        case Failure(_: UnsupportedOperationException) => //we're not testing tsv compatibility
        case Failure(t) => throw t
      })
    }

    //construct a test for both over and under
    List(
      s"read $command" should "limit according to a setting" in testOver,
      it should "allow when under the  limit" in testUnder
    )
  }

  //test all the functions
  List[(String, Int, ReadLikeFunctions => (Seq[Try[WdlValue]] => Try[WdlValue]))](
    ("lines", _readLinesLimit, _.read_lines),
    ("int", _readIntLimit, _.read_int),
    ("map", _readMapLimit, _.read_map),
    ("float", _readFloatLimit, _.read_float),
    ("String", _readStringLimit, _.read_string),
    ("tsv", _readTsvLimit, _.read_tsv),
    ("object", _readObjectLimit, _.read_object)
  ).flatMap {
    (testOverUnder _).tupled
  }
}
