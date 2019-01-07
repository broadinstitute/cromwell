package centaur

import centaur.test.standard.CentaurTestCase
import org.scalatest.DoNotDiscover

@DoNotDiscover
class PapiUpgradeTestCaseSpec(cromwellBackends: List[String])
  extends AbstractCromwellEngineOrBackendUpgradeTestCaseSpec(cromwellBackends) {

  def this() = this(CentaurTestSuite.cromwellBackends)

  override def testType: String = "PAPI upgrade"

  override def isMatchingUpgradeTest(testCase: CentaurTestCase): Boolean = CentaurTestSuite.isPapiUpgradeTest(testCase)
}
