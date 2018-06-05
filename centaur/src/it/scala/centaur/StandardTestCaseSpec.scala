package centaur

import org.scalatest._

@DoNotDiscover
class StandardTestCaseSpec(cromwellBackends: List[String]) extends AbstractCentaurTestCaseSpec(cromwellBackends) with ParallelTestExecution {
  
  def this() = this(CentaurTestSuite.cromwellBackends)

  allTestCases.filter(CentaurTestSuite.runParallel) foreach executeStandardTest

  // The WDL version upgrade tests are just regular draft-2 test cases tagged for re-use in testing the upgrade script
  allTestCases.filter(CentaurTestSuite.runUpgradeTest) foreach executeUpgradeTest

}
