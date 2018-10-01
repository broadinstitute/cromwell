package centaur

import org.scalatest.{DoNotDiscover, ParallelTestExecution}

@DoNotDiscover
class WdlUpgradeTestCaseSpec(cromwellBackends: List[String])
  extends AbstractCentaurTestCaseSpec(cromwellBackends) with ParallelTestExecution with CentaurTestSuiteShutdown {

  def this() = this(CentaurTestSuite.cromwellBackends)

  // The WDL version upgrade tests are just regular draft-2 test cases tagged for re-use in testing the upgrade script
  allTestCases.filter(CentaurTestSuite.isWdlUpgradeTest) foreach executeWdlUpgradeTest
}
