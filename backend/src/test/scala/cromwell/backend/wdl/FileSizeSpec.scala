package cromwell.backend.wdl

import java.nio.file.StandardOpenOption._
import java.nio.file.{Path, Paths}

import cats.effect.IO
import com.google.common.io.Files
import cromwell.backend.standard.{DefaultStandardExpressionFunctionsParams, StandardExpressionFunctions}
import cromwell.core.Tags.PostWomTest
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.{CallContext, TestKitSuite}
import fs2.Stream
import org.scalatest.{FlatSpecLike, Matchers}
import wom.values._

import scala.util.{Failure, Success, Try}

class FileSizeSpec extends TestKitSuite with FlatSpecLike with Matchers {
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

    val dp = DefaultStandardExpressionFunctionsParams(
      List(cromwell.core.path.DefaultPathBuilder),
      CallContext(path, "stdout", "stderr"),
      simpleIoActor,
      scala.concurrent.ExecutionContext.global)

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

  def testOverUnder(command: String, n: Int, f: ReadLikeFunctions => (Seq[Try[WomValue]] => Try[WomValue])) = {

    def testInner(n: Int, test: Try[WomValue] => Unit) = {

      def createTempFileOfSize(size: Int): Path = {

        val fn = tempDir.toString + "/" + scala.util.Random.alphanumeric.take(5).mkString
        val jPath = Paths.get(fn)
        jPath.toFile.deleteOnExit
        val start = Stream.eval[IO, Byte](IO(1.toByte)).repeat.take(size.toLong)
        val end = fs2.io.file.writeAll[IO](jPath, Seq(CREATE_NEW, WRITE))
        (start to end).run.unsafeRunSync
        //jPath is now a file of n bytes, we can return it
        jPath
      }

      val file = createTempFileOfSize(n)
      val params = Seq(Try(WomString(file.toString)))

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
      s"read $command" should "limit according to a setting" taggedAs PostWomTest ignore testOver,
      it should "allow when under the  limit" in testUnder
    )
  }

  //test all the functions
  // TODO WOM: Restore those tests. Should they be in wdl4s ? Hard to tell without https://github.com/broadinstitute/cromwell/issues/2611
  List[(String, Int, ReadLikeFunctions => (Seq[Try[WomValue]] => Try[WomValue]))](
//    ("lines", _readLinesLimit, _.read_lines),
//    ("int", _readIntLimit, _.read_int),
//    ("map", _readMapLimit, _.read_map),
//    ("float", _readFloatLimit, _.read_float),
//    ("String", _readStringLimit, _.read_string),
//    ("tsv", _readTsvLimit, _.read_tsv),
//    ("object", _readObjectLimit, _.read_object)
  ).flatMap {
    (testOverUnder _).tupled
  }
}
