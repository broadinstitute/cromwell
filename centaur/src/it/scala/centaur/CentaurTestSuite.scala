package centaur

import centaur.api.CentaurCromwellClient
import centaur.test.standard.CentaurTestCase
import centaur.test.standard.CentaurTestFormat.{InstantAbort, RestartFormat, ScheduledAbort}
import org.scalatest.{BeforeAndAfterAll, ParallelTestExecution, Suite, Suites}

import scala.sys.ShutdownHookThread

object CentaurTestSuite {
  // Start cromwell if we're in Managed mode
  // Note: we can't use beforeAll to start Cromwell, because beforeAll is executed once the suite is instantiated and the
  // tests exist. However because the set of tests differs depending on the backends supported by Cromwell, it needs to be up
  // before we can generate the tests.
  startCromwell()

  def startCromwell(): Unit = {
    CentaurConfig.runMode match {
      case ManagedCromwellServer(preRestart, _, _) =>
        CromwellManager.startCromwell(preRestart)
      case _ =>
    }
  }

  val cromwellBackends = CentaurCromwellClient.backends.get.supportedBackends.map(_.toLowerCase)
  
  def runSequential(testCase: CentaurTestCase) = testCase.testFormat match {
    case _: RestartFormat| _: ScheduledAbort | InstantAbort => true
    case _ => false
  }

  def isWdlUpgradeTest(testCase: CentaurTestCase): Boolean = testCase.containsTag("wdl_upgrade")

  def isEngineUpgradeTest(testCase: CentaurTestCase): Boolean = testCase.containsTag("engine_upgrade")

  def isPapiUpgradeTest(testCase: CentaurTestCase): Boolean = testCase.containsTag("papi_upgrade")

  def runParallel(testCase: CentaurTestCase) = !runSequential(testCase)
}

/**
  * Shuts down the managed cromwell after centaur finishes running.
  * Should be mixed into any root suite of centaur tests.
  */
trait CentaurTestSuiteShutdown extends Suite with BeforeAndAfterAll {
  private var shutdownHook: Option[ShutdownHookThread] = _

  override protected def beforeAll() = {
    shutdownHook = Option(sys.addShutdownHook { CromwellManager.stopCromwell("JVM Shutdown Hook") })
  }

  override protected def afterAll() = {
    CromwellManager.stopCromwell("ScalaTest AfterAll")
    shutdownHook.foreach(_.remove())
  }
}

/**
  * The main centaur test suites, runs sub suites in parallel, but allows better control over the way each nested suite runs.
  */
class CentaurTestSuite
  extends Suites(new SequentialTestCaseSpec(), new StandardTestCaseSpec())
    with ParallelTestExecution
    with CentaurTestSuiteShutdown
