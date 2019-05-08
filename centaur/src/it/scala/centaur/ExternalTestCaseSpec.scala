package centaur

import better.files._
import cats.data.Validated.{Invalid, Valid}
import centaur.test.standard.CentaurTestCase
import com.typesafe.scalalogging.StrictLogging

class ExternalTestCaseSpec(cromwellBackends: List[String]) extends AbstractCentaurTestCaseSpec(cromwellBackends) with StrictLogging {

  def this() = this(CentaurTestSuite.cromwellBackends)

  val testFile = sys.env.get("CENTAUR_TEST_FILE") match {
    case Some(filePath) => runTestFile(filePath)
    case _ =>
      logger.info("No external test to run")
  }

  def runTestFile(testFile: String) = {
    CentaurTestCase.fromFile(cromwellTracker = None)(File(testFile)) match {
      case Valid(testCase) => executeStandardTest(testCase)
      case Invalid(error) =>
        fail(s"Invalid test case: ${error.toList.mkString(", ")}")
    }
  }
}
