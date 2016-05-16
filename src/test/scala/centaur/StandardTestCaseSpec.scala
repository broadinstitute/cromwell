package centaur

import java.nio.file.Path

import cats.data.Validated.{Invalid, Valid}
import cats.syntax.traverse._
import cats.std.list._
import centaur.test.Test
import centaur.test.standard.StandardTestCase

import scala.language.postfixOps
import org.scalatest._

class StandardTestCaseSpec extends FlatSpec with Matchers with ParallelTestExecution{
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
    case t => executeStandardTest(t, t.testFormat.testFunction)
  }

  def executeStandardTest(testCase: StandardTestCase, f: WorkflowRequest => Test[_]): Unit = {
    it should s"${testCase.testFormat.testSpecString} ${testCase.name}" in {
      f(WorkflowRequest(testCase)).run.get
    }
  }
}
