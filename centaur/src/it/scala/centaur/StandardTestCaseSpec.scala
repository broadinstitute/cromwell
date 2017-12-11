package centaur

import org.scalatest._

@DoNotDiscover
class StandardTestCaseSpec(cromwellBackends: List[String]) extends AbstractCentaurTestCaseSpec(cromwellBackends) with ParallelTestExecution {
  
  def this() = this(CentaurTestSuite.cromwellBackends)

  allTestCases.filter(CentaurTestSuite.runParallel) foreach executeStandardTest
}
