package centaur

import java.nio.file.Path

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.traverse._
import cats.std.list._
import centaur.test.Test
import centaur.test.standard.StandardTestCase
import centaur.test.workflow.Workflow

import scala.language.postfixOps
import org.scalatest._

class StandardTestCaseSpec extends FlatSpec with Matchers with ParallelTestExecution {
  def testCases(basePath: Path): List[StandardTestCase] = {
    // IntelliJ will give some red squiggles in the following block. It lies.
    basePath.toFile.listFiles.toList collect { case x if x.isFile => x.toPath } traverse StandardTestCase.fromPath match {
      case Valid(l) => l
      case Invalid(e) => throw new IllegalStateException("\n" + e.unwrap.mkString("\n") + "\n")
    }
  }

  // Optional test cases are provided by the end user as opposed to the ones built in to the system
  val optionalTestCases = CentaurConfig.optionalTestPath map testCases getOrElse List.empty
  optionalTestCases ++ testCases(CentaurConfig.standardTestCasePath) foreach {
    case t => executeStandardTest(t, t.testFunction)
  }

  def executeStandardTest(testCase: StandardTestCase, f: Workflow => Test[_]): Unit = {

    def nameTest = it should s"${testCase.testFormat.testSpecString} ${testCase.workflow.name}"
    def runTest = f(testCase.workflow).run.get

    // Make tags, but enforce lowercase:
    val tags = (testCase.tagStrings :+ testCase.workflow.name :+ testCase.testFormat.name) map { x => Tag(x.toLowerCase) }

    tags match {
      case Nil => nameTest in runTest
      case head :: Nil => nameTest taggedAs head in runTest
      case head :: tail => nameTest taggedAs(head, tail: _*) in runTest
    }
  }
}
