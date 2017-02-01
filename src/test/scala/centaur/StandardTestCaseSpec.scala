package centaur

import java.nio.file.Path

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.traverse._
import centaur.api.CromwellBackendsCompanion
import centaur.test.{ErrorOr, Test}
import centaur.test.standard.StandardTestCase
import centaur.test.workflow.Workflow

import scala.language.postfixOps
import org.scalatest._

class StandardTestCaseSpec extends FlatSpec with Matchers with ParallelTestExecution {
  def testCases(basePath: Path): List[StandardTestCase] = {
    val files = basePath.toFile.listFiles.toList collect { case x if x.isFile => x.toPath }
    val testCases = files.traverse[ErrorOr, StandardTestCase](StandardTestCase.fromPath)

    testCases match {
      case Valid(l) => l
      case Invalid(e) => throw new IllegalStateException("\n" + e.toList.mkString("\n") + "\n")
    }
  }

  val cromwellBackends = CromwellBackendsCompanion.supportedBackends

  // Optional test cases are provided by the end user as opposed to the ones built in to the system
  val optionalTestCases = CentaurConfig.optionalTestPath map testCases getOrElse List.empty
  val allTestCases = optionalTestCases ++ testCases(CentaurConfig.standardTestCasePath)
  allTestCases foreach { t => executeStandardTest(t, t.testFunction) }

  def executeStandardTest(testCase: StandardTestCase, f: Workflow => Test[_]): Unit = {
    def nameTest = it should s"${testCase.testFormat.testSpecString} ${testCase.workflow.testName}"
    def runTest = f(testCase.workflow).run.get

    // Make tags, but enforce lowercase:
    val tags = (testCase.testOptions.tags :+ testCase.workflow.testName :+ testCase.testFormat.name) map { x => Tag(x.toLowerCase) }
    val isIgnored = testCase.isIgnored(cromwellBackends)

    tags match {
      case Nil => runOrDont(nameTest, isIgnored, runTest)
      case head :: Nil => runOrDont(nameTest taggedAs head, isIgnored, runTest)
      case head :: tail => runOrDont(nameTest taggedAs(head, tail: _*), isIgnored, runTest)
    }
  }

  private def runOrDont(itVerbString: ItVerbString, ignore: Boolean, runTest: => Any) = {
    if (ignore) {
      itVerbString ignore runTest
    } else {
      itVerbString in runTest
    }
  }

  private def runOrDont(itVerbStringTaggedAs: ItVerbStringTaggedAs, ignore: Boolean, runTest: => Any) = {
    if (ignore) {
      itVerbStringTaggedAs ignore runTest
    } else {
      itVerbStringTaggedAs in runTest
    }
  }
}
