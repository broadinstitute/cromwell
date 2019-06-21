package centaur

import centaur.api.CentaurCromwellClient
import centaur.test.standard.CentaurTestCase
import com.typesafe.config.{Config, ConfigFactory}
import com.typesafe.scalalogging.StrictLogging
import net.ceedubs.ficus.Ficus._
import org.scalatest.{BeforeAndAfterAll, ParallelTestExecution, Suite, Suites}

import scala.sys.ShutdownHookThread

object CentaurTestSuite extends StrictLogging {
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

  val cromwellBackends = CentaurCromwellClient.backends.unsafeRunSync().supportedBackends.map(_.toLowerCase)

  def isWdlUpgradeTest(testCase: CentaurTestCase): Boolean = testCase.containsTag("wdl_upgrade")

  def isEngineUpgradeTest(testCase: CentaurTestCase): Boolean = testCase.containsTag("engine_upgrade")

  def isPapiUpgradeTest(testCase: CentaurTestCase): Boolean = testCase.containsTag("papi_upgrade")

  /** Horicromtality-related assertion config. */
  val cromwellTracker: Option[CromwellTracker] = {

    def backendCountFromConfig(config: Config): Option[Int] = {
      val assert = config.getOrElse("assert", default = false)
      val backendCount = config.as[Option[Int]]("backend-count")
      (assert, backendCount) match {
        case (false, _) => None
        case (true, Some(_)) => backendCount
        case (true, _) =>
          val message = "Invalid Centaur configuration: `horicromtal` must define `backend-count` if `assert = true`"
          throw new RuntimeException(message)
      }
    }

    for {
      config <- ConfigFactory.load().as[Option[Config]]("centaur.horicromtal")
      backendCount <- backendCountFromConfig(config)
      configuredSignificance = config.getOrElse("significance-level", 0.05)
    } yield CromwellTracker(backendCount, configuredSignificance)
  }
  logger.info(s"Horicromtal tracker config: {}", cromwellTracker)
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
    CentaurTestSuite.cromwellTracker foreach { _.assertHoricromtality() }
    shutdownHook.foreach(_.remove())
  }
}

/**
  * The main centaur test suites, runs sub suites in parallel, but allows better control over the way each nested suite runs.
  */
class CentaurTestSuite
  extends Suites(new SequentialTestCaseSpec(), new ParallelTestCaseSpec())
    with ParallelTestExecution
    with CentaurTestSuiteShutdown
