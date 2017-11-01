package centaur

import java.nio.file.Path

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.traverse._
import centaur.test.standard.CentaurTestCase
import lenthall.validation.ErrorOr.ErrorOr
import org.scalatest.{DoNotDiscover, FlatSpec, Matchers, Tag}

@DoNotDiscover
abstract class AbstractCentaurTestCaseSpec(cromwellBackends: List[String]) extends FlatSpec with Matchers {

  private def testCases(basePath: Path): List[CentaurTestCase] = {
    val files = basePath.toFile.listFiles.toList collect { case x if x.isFile => x.toPath }
    val testCases = files.traverse[ErrorOr, CentaurTestCase](CentaurTestCase.fromPath)

    testCases match {
      case Valid(l) => l
      case Invalid(e) => throw new IllegalStateException("\n" + e.toList.mkString("\n") + "\n")
    }
  }
  
  def allTestCases: List[CentaurTestCase] = {
    val optionalTestCases = CentaurConfig.optionalTestPath map testCases getOrElse List.empty
    val standardTestCases = testCases(CentaurConfig.standardTestCasePath)
    optionalTestCases ++ standardTestCases
  }

  def executeStandardTest(testCase: CentaurTestCase): Unit = {
    def nameTest = s"${testCase.testFormat.testSpecString} ${testCase.workflow.testName}"
    def runTest = testCase.testFunction.run.get

    // Make tags, but enforce lowercase:
    val tags = (testCase.testOptions.tags :+ testCase.workflow.testName :+ testCase.testFormat.name) map { x => Tag(x.toLowerCase) }
    val isIgnored = testCase.isIgnored(cromwellBackends)

    runOrDont(nameTest, tags, isIgnored, runTest)
  }

  private def runOrDont(testName: String, tags: List[Tag], ignore: Boolean, runTest: => Any): Unit = {

    val itShould: ItVerbString = it should testName

    tags match {
      case Nil => runOrDont(itShould, ignore, runTest)
      case head :: Nil => runOrDont(itShould taggedAs head, ignore, runTest)
      case head :: tail => runOrDont(itShould taggedAs(head, tail: _*), ignore, runTest)
    }
  }

  private def runOrDont(itVerbString: ItVerbString, ignore: Boolean, runTest: => Any): Unit = {
    if (ignore) {
      itVerbString ignore runTest
    } else {
      itVerbString in runTest
    }
  }

  private def runOrDont(itVerbStringTaggedAs: ItVerbStringTaggedAs, ignore: Boolean, runTest: => Any): Unit = {
    if (ignore) {
      itVerbStringTaggedAs ignore runTest
    } else {
      itVerbStringTaggedAs in runTest
    }
  }
  
}
