import Dependencies._
import sbt.Defaults._
import sbt.Keys._
import sbt._
import complete.DefaultParsers._
import sbt.util.Logger

import scala.sys.process._

object Testing {
  private val AllTests = config("alltests") extend Test
  private val CromwellBenchmarkTest = config("benchmark") extend Test

  private val DockerTestTag = "DockerTest"
  private val CromwellIntegrationTestTag = "CromwellIntegrationTest"
  private val GcsIntegrationTestTag = "GcsIntegrationTest"
  private val AwsIntegrationTestTag = "AwsTest"
  private val DbmsTestTag = "DbmsTest"

  private val AllTestTags = List(
    DockerTestTag,
    CromwellIntegrationTestTag,
    GcsIntegrationTestTag,
    AwsIntegrationTestTag,
    DbmsTestTag
  )

  val minnieKenny = inputKey[Unit]("Run minnie-kenny.")

  private val excludeTestTags: Seq[String] =
    sys.env
      .get("CROMWELL_SBT_TEST_EXCLUDE_TAGS")
      .filter(_.nonEmpty)
      .map(_.split(",").toList.map(_.trim))
      .getOrElse(AllTestTags)

  private val spanScaleFactor: String = sys.env.getOrElse("CROMWELL_SBT_TEST_SPAN_SCALE_FACTOR", "1")

  /*
  The arguments that will be added to the default test config, but removed from all other configs.
  `sbt coverage test` adds other arguments added to generate the coverage reports.
  Tracking the arguments we add to the default allows one to later remove them when building up other configurations.
 */
  private val excludeTestArgs = excludeTestTags.map(Tests.Argument(TestFrameworks.ScalaTest, "-l", _))

  private val TestReportArgs =
    Tests.Argument(
      TestFrameworks.ScalaTest,
      "-oDSI",
      "-h",
      "target/test-reports",
      "-u",
      "target/test-reports",
      "-F",
      spanScaleFactor,
      "-W",
      "300",
      "300",
    )

  /** Run minnie-kenny only once per sbt invocation. */
  class MinnieKennySingleRunner() {
    private val mutex = new Object
    private var resultOption: Option[Int] = None

    /** Run using the logger, throwing an exception only on the first failure. */
    def runOnce(log: Logger, args: Seq[String]): Unit = {
      mutex synchronized {
        if (resultOption.isEmpty) {
          log.debug(s"Running minnie-kenny.sh${args.mkString(" ", " ", "")}")
          val result = ("./minnie-kenny.sh" +: args) ! log
          resultOption = Option(result)
          if (result == 0)
            log.debug("Successfully ran minnie-kenny.sh")
          else
            sys.error("Running minnie-kenny.sh failed. Please double check for errors above.")
        }
      }
    }
  }

  // Only run one minnie-kenny.sh at a time!
  private lazy val minnieKennySingleRunner = new MinnieKennySingleRunner

  val testSettings = List(
    libraryDependencies ++= testDependencies.map(_ % Test),
    // `test` (or `assembly`) - Run all tests, except docker and integration and DBMS
    testOptions in Test ++= Seq(TestReportArgs) ++ excludeTestArgs,
    // `alltests:test` - Run all tests
    testOptions in AllTests := (testOptions in Test).value.diff(excludeTestArgs),
    // Add scalameter as a test framework in the CromwellBenchmarkTest scope
    testFrameworks in CromwellBenchmarkTest += new TestFramework("org.scalameter.ScalaMeterFramework"),
    // Don't execute benchmarks in parallel
    parallelExecution in CromwellBenchmarkTest := false,
    // Make sure no secrets are commited to git
    minnieKenny := {
      val log = streams.value.log
      val args = spaceDelimited("<arg>").parsed
      minnieKennySingleRunner.runOnce(log, args)
    },
    test in Test := {
      minnieKenny.toTask("").value
      (test in Test).value
    },
  )

  val integrationTestSettings = List(
    libraryDependencies ++= testDependencies.map(_ % IntegrationTest)
  ) ++ itSettings

  def addTestSettings(project: Project) = {
    project
      .settings(testSettings)
      .configs(AllTests).settings(inConfig(AllTests)(Defaults.testTasks): _*)
      .configs(CromwellBenchmarkTest).settings(inConfig(CromwellBenchmarkTest)(Defaults.testTasks): _*)
  }

  def addIntegrationTestSettings(project: Project) = {
    project
      .settings(integrationTestSettings)
      .configs(IntegrationTest)
  }

}
