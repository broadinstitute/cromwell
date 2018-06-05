package centaur

import org.scalatest.ParallelTestExecution

class UpgradeTestCaseSpec(cromwellBackends: List[String]) extends AbstractCentaurTestCaseSpec(cromwellBackends) with ParallelTestExecution {

  def this() = this(CentaurTestSuite.cromwellBackends)

  allTestCases.filter(CentaurTestSuite.runUpgradeTest) foreach executeUpgradeTest

}
