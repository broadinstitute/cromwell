package centaur

import java.nio.file.Path

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.traverse._
import centaur.test.standard.CentaurTestCase
import cromwell.core.path
import cromwell.core.path.DefaultPathBuilder
import org.scalatest._

import scala.concurrent.Future

@DoNotDiscover
abstract class AbstractCentaurTestCaseSpec(cromwellBackends: List[String]) extends AsyncFlatSpec with Matchers {

  private def testCases(basePath: Path): List[CentaurTestCase] = {
    val files = basePath.toFile.listFiles.toList collect { case x if x.isFile => x.toPath }
    val testCases = files.traverse(CentaurTestCase.fromPath)

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
    def runTest() = {
      testCase.testFunction.run.unsafeToFuture().map(_ => assert(true))
    }

    // Make tags, but enforce lowercase:
    val tags = (testCase.testOptions.tags :+ testCase.workflow.testName :+ testCase.testFormat.name) map { x => Tag(x.toLowerCase) }
    val isIgnored = testCase.isIgnored(cromwellBackends)

    runOrDont(nameTest, tags, isIgnored, runTest())
  }

  def executeUpgradeTest(testCase: CentaurTestCase): Unit =
    executeStandardTest(upgradeTestWdl(testCase))

  private def upgradeTestWdl(testCase: CentaurTestCase): CentaurTestCase = {
    import womtool.WomtoolMain

    // The suffix matters because WomGraphMaker.getBundle() uses it to choose the language factory
    val tempFile: path.Path = DefaultPathBuilder.createTempFile(suffix = "wdl").append(testCase.workflow.data.workflowContent.get) //TODO: Saloni-What about this?
    val upgradeResult = WomtoolMain.upgrade(tempFile.pathAsString)

    testCase.copy(
      workflow = testCase.workflow.copy(
        testName = testCase.workflow.testName + " (draft-2 to 1.0 upgrade)",
        data = testCase.workflow.data.copy(
          workflowContent = Option(upgradeResult.stdout.get)))) ////TODO: Saloni-What about this?
  }

  private def runOrDont(testName: String, tags: List[Tag], ignore: Boolean, runTest: => Future[Assertion]): Unit = {

    val itShould: ItVerbString = it should testName

    tags match {
      case Nil => runOrDont(itShould, ignore, runTest)
      case head :: Nil => runOrDont(itShould taggedAs head, ignore, runTest)
      case head :: tail => runOrDont(itShould taggedAs(head, tail: _*), ignore, runTest)
    }
  }

  private def runOrDont(itVerbString: ItVerbString, ignore: Boolean, runTest: => Future[Assertion]): Unit = {
    if (ignore) {
      itVerbString ignore runTest
    } else {
      itVerbString in runTest
    }
  }

  private def runOrDont(itVerbStringTaggedAs: ItVerbStringTaggedAs, ignore: Boolean, runTest: => Future[Assertion]): Unit = {
    if (ignore) {
      itVerbStringTaggedAs ignore runTest
    } else {
      itVerbStringTaggedAs in runTest
    }
  }
  
}
