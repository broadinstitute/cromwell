package centaur

import org.scalatest.{DoNotDiscover, ParallelTestExecution}

@DoNotDiscover
class UpgradeTestCaseSpec(cromwellBackends: List[String]) extends AbstractCentaurTestCaseSpec(cromwellBackends) with ParallelTestExecution {

  def this() = this(CentaurTestSuite.cromwellBackends)

  // The WDL version upgrade tests are just regular draft-2 test cases tagged for re-use in testing the upgrade script
  allTestCases.filter(CentaurTestSuite.isUpgradeTest) foreach executeUpgradeTest
}
