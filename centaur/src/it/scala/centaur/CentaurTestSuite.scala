package centaur

import centaur.api.CentaurCromwellClient
import centaur.test.standard.CentaurTestCase
import centaur.test.standard.CentaurTestFormat.{InstantAbort, RestartFormat, ScheduledAbort}
import org.scalatest.{BeforeAndAfterAll, ParallelTestExecution, Suites}

import scala.sys.ShutdownHookThread

object CentaurTestSuite {
  // Start cromwell if we're in Managed mode
  // Note: we can't use beforeAll to start Cromwell, because beforeAll is executed once the suite is instantiated and the
  // tests exist. However because the set of tests differs depending on the backends supported by Cromwell, it needs to be up
  // before we can generate the tests.
  CentaurConfig.runMode match {
    case ManagedCromwellServer(preRestart, _, _) => 
      CromwellManager.startCromwell(preRestart)
    case _ =>
  }

  val cromwellBackends = CentaurCromwellClient.backends.get.supportedBackends.map(_.toLowerCase)
  
  def runSequential(testCase: CentaurTestCase) = testCase.testFormat match {
    case _: RestartFormat| _: ScheduledAbort | InstantAbort => true
    case _ => false
  }
  
  def runParallel(testCase: CentaurTestCase) = !runSequential(testCase)
}

/**
  * The main centaur test suites, runs sub suites in parallel, but allows better control over the way each nested suite runs.
  */
class CentaurTestSuite extends Suites(new SequentialTestCaseSpec(), new StandardTestCaseSpec()) with ParallelTestExecution with BeforeAndAfterAll {
  private var shutdownHook: Option[ShutdownHookThread] = _

  override def beforeAll() = {
    shutdownHook = Option(sys.addShutdownHook { CromwellManager.stopCromwell() })
  }

  override def afterAll() = {
    CromwellManager.stopCromwell()
    shutdownHook.foreach(_.remove())
  }
}
