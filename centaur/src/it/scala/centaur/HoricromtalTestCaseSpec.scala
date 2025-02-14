package centaur

import centaur.test.standard.CentaurTestCase
import org.scalatest.DoNotDiscover

@DoNotDiscover
class HoricromtalTestCaseSpec(cromwellBackends: List[String]) extends AbstractCentaurTestCaseSpec(cromwellBackends) {

  def this() = this(CentaurTestSuite.cromwellBackends)

  def testType: String = "Horicromtal"

  def isMatchingHoricromtalTest(testCase: CentaurTestCase): Boolean = CentaurTestSuite.isHoricromtalTest(testCase)
}
