package centaur

import org.scalatest.{DoNotDiscover, Matchers}

/**
  * All test cases that trigger a Cromwell restart. Note that this suite does not mix in ParallelTestExecution
  * such that the restarting tests execute sequentially to avoid a mayhem of Cromwell restarts
 */
@DoNotDiscover
class SequentialTestCaseSpec(cromwellBackends: List[String]) extends AbstractCentaurTestCaseSpec(cromwellBackends) with Matchers {

  def this() = this(CentaurTestSuite.cromwellBackends)

  allTestCases.filterNot(_.testFormat.isParallel) foreach executeStandardTest

}
