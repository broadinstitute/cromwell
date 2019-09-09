import Dependencies._
import sbt.Defaults._
import sbt.Keys._
import sbt._

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

  val testSettings = List(
    libraryDependencies ++= testDependencies.map(_ % Test),
    // `test` (or `assembly`) - Run all tests, except docker and integration and DBMS
    testOptions in Test ++= Seq(TestReportArgs) ++ excludeTestArgs,
    // `alltests:test` - Run all tests
    testOptions in AllTests := (testOptions in Test).value.diff(excludeTestArgs),
    // Add scalameter as a test framework in the CromwellBenchmarkTest scope
    testFrameworks in CromwellBenchmarkTest += new TestFramework("org.scalameter.ScalaMeterFramework"),
    // Don't execute benchmarks in parallel
    parallelExecution in CromwellBenchmarkTest := false
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
