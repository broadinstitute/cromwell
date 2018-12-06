package centaur

import centaur.test.standard.CentaurTestCase
import org.scalatest.DoNotDiscover

@DoNotDiscover
class EngineUpgradeTestCaseSpec(cromwellBackends: List[String]) extends
  AbstractCromwellEngineOrBackendUpgradeTestCaseSpec(cromwellBackends) {

  def this() = this(CentaurTestSuite.cromwellBackends)

  override def testType: String = "Engine upgrade"

  override def isMatchingUpgradeTest(testCase: CentaurTestCase): Boolean = CentaurTestSuite.isEngineUpgradeTest(testCase)
}
