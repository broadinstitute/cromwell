package cromwell.backend.wdl

import java.nio.file.{Paths, Path}
import java.nio.file.StandardOpenOption._
import scala.util.{Failure, Success, Try}

import cromwell.backend.standard.{DefaultStandardExpressionFunctionsParams, StandardExpressionFunctions}
import cromwell.core.CallContext
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.values._
import com.google.common.io.Files
import fs2.{Task, Stream}
import java.nio.file.Files.delete

class FileSizeSpec extends FlatSpec with Matchers /*with TryValues*/ {
  final val _readLinesLimit = 4
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

    val dp = DefaultStandardExpressionFunctionsParams(List(cromwell.core.path.DefaultPathBuilder), CallContext(path, "stdout", "stderr" ))

    new StandardExpressionFunctions(dp) {
      override val readLinesLimit = _readLinesLimit
      override val readIntLimit = _readIntLimit
      override val readFloatLimit = _readFloatLimit
      override val readStringLimit = _readStringLimit
      override val readJsonLimit = _readJsonLimit
      override val readTsvLimit = _readTsvLimit
      override val readMapLimit = _readMapLimit
      override val readObjectLimit = _readObjectLimit
    }
  }

  def createTempFileOfSize(n: Int): Path = {
    val tempDir = Files.createTempDir

    val fn = tempDir.toString + "/whatever"

    val jPath = Paths.get(fn)

    val end = fs2.io.file.writeAll[Task](jPath, Seq(CREATE_NEW, WRITE) )

    val start = Stream[Task, Byte](1).repeat.take(n.toLong)

    (start to end).run.unsafeRunSync 
    jPath
  }

  def testInner(n: Int, f: ReadLikeFunctions => (Seq[Try[WdlValue]] =>Try[WdlValue]),g: Try[WdlValue] => Unit) = {

    val file = createTempFileOfSize(n)

    val params = Seq(Try(WdlString(file.toString)))

    f(rlf)(params) match {
      case t => g(t)
    }

    delete(file)
  }

  def testOver(n: Int, f: ReadLikeFunctions => (Seq[Try[WdlValue]] =>Try[WdlValue])) = {
    testInner(n + 1, f, {
      case Failure(s: FileSizeTooBig) => //success
      case t => throw new RuntimeException(s"should not have eaten this file that is too big! msg: $t")
    })
  }

  def testUnder(n: Int, f: ReadLikeFunctions => (Seq[Try[WdlValue]] =>Try[WdlValue])) = {
    testInner(n - 1, f, {
      case Success(_) => 
      case Failure(t) => throw t
    })
  }

  "read lines" should "limit based upon argument" in 
    testOver(_readLinesLimit + 1, _.read_lines)

  "read lines" should "allow based upon argument" in 
    testUnder(_readLinesLimit - 1, _.read_lines)

  "read int" should "limit based upon argument" in 
    testOver(_readIntLimit + 1, _.read_int)

    /*
  "read lines" should "allow based upon argument" in 
    testUnder(_readLinesLimit - 1, _.read_lines)
    */
}
