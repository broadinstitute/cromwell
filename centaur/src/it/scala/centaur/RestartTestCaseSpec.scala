package centaur

import centaur.test.standard.CentaurTestFormat.RestartFormat
import org.scalatest.{DoNotDiscover, Matchers}

/**
  * All test cases that trigger a Cromwell restart. Note that this suite does not mix in ParallelTestExecution
  * such that the restarting tests execute sequentially to avoid a mayhem of Cromwell restarts
 */
@DoNotDiscover
class RestartTestCaseSpec(cromwellBackends: List[String]) extends AbstractCentaurTestCaseSpec(cromwellBackends) with Matchers {

  def this() = this(CentaurTestSuite.cromwellBackends)
  
  allTestCases.filter( _.testFormat match {
    case _: RestartFormat => true
    case _ => false
  }
  ) foreach executeStandardTest

}
