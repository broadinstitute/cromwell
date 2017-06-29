package centaur

import java.nio.file.Path

import cats.data.Validated.{Invalid, Valid}
import cats.instances.list._
import cats.syntax.traverse._
import centaur.api.CentaurCromwellClient
import centaur.test.ErrorOr
import centaur.test.standard.CentaurTestCase

import scala.language.postfixOps
import org.scalatest._

class StandardTestCaseSpec extends FlatSpec with Matchers with ParallelTestExecution with BeforeAndAfterAll {
  
  // Start cromwell if we're in Managed mode
  // Note: we can't use beforeAll to start Cromwell, because beforeAll is executed once the suite is instantiated and the
  // tests exist. However because the set of tests differs depending on the backends supported by Cromwell, it needs to be up
  // before we can generate the tests.
  // The solution chosen is to use a singleton object containing the Cromwell server state and ways to start / stop it.
  // Another possibly better way would be to use something like https://stackoverflow.com/a/15556379/1498572 ?
  CentaurConfig.runMode match {
    case ManagedCromwellServer(preRestart, _, _) => CromwellManager.startCromwell(preRestart)
    case _ =>
  }
  
  private val cromwellBackends = CentaurCromwellClient.backends.get.supportedBackends.map(_.toLowerCase)

  override def beforeAll() = {
      sys.addShutdownHook { CromwellManager.stopCromwell() }
  }
  
  override def afterAll() = CromwellManager.stopCromwell()

  def testCases(basePath: Path): List[CentaurTestCase] = {
    val files = basePath.toFile.listFiles.toList collect { case x if x.isFile => x.toPath }
    val testCases = files.traverse[ErrorOr, CentaurTestCase](CentaurTestCase.fromPath)

    testCases match {
      case Valid(l) => l
      case Invalid(e) => throw new IllegalStateException("\n" + e.toList.mkString("\n") + "\n")
    }
  }

  // Optional test cases are provided by the end user as opposed to the ones built in to the system
  private val optionalTestCases = CentaurConfig.optionalTestPath map testCases getOrElse List.empty
  private val standardTestCases = testCases(CentaurConfig.standardTestCasePath)
  (optionalTestCases ++ standardTestCases) foreach executeStandardTest

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
