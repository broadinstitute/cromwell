package centaur

import centaur.CentaurTestSuite.cromwellTracker
import com.typesafe.scalalogging.LazyLogging
import org.scalatest._


/**
  * Runs test cases in parallel, this should be the default type for tests unless they would otherwise crosstalk in undesirable
  * ways with other tests and must be made sequential.
  */
@DoNotDiscover
class ParallelTestCaseSpec(cromwellBackends: List[String])
  extends AbstractCentaurTestCaseSpec(cromwellBackends, cromwellTracker = cromwellTracker) with ParallelTestExecution with LazyLogging {
  
  def this() = this(CentaurTestSuite.cromwellBackends)

  allTestCases.filter(_.testFormat.isParallel) foreach { test =>
    import java.time
    logger.info(s"Starting test ${test.name} at ${time.LocalDateTime.now().toString}")
    executeStandardTest(test)
  }
}
