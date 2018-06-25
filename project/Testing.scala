import Dependencies._
import sbt.Defaults._
import sbt.Keys._
import sbt._

object Testing {
  lazy val AllTests = config("alltests") extend Test
  lazy val DockerTest = config("docker") extend Test
  lazy val NoDockerTest = config("nodocker") extend Test
  lazy val CromwellIntegrationTest = config("integration") extend Test
  lazy val CromwellBenchmarkTest = config("benchmark") extend Test
  lazy val CromwellNoIntegrationTest = config("nointegration") extend Test
  lazy val DbmsTest = config("dbms") extend Test

  lazy val DockerTestTag = "DockerTest"
  lazy val UseDockerTaggedTests = Tests.Argument(TestFrameworks.ScalaTest, "-n", DockerTestTag)
  lazy val DontUseDockerTaggedTests = Tests.Argument(TestFrameworks.ScalaTest, "-l", DockerTestTag)

  lazy val CromwellIntegrationTestTag = "CromwellIntegrationTest"
  lazy val UseCromwellIntegrationTaggedTests = Tests.Argument(TestFrameworks.ScalaTest, "-n", CromwellIntegrationTestTag)
  lazy val DontUseCromwellIntegrationTaggedTests = Tests.Argument(TestFrameworks.ScalaTest, "-l", CromwellIntegrationTestTag)

  lazy val GcsIntegrationTestTag = "GcsIntegrationTest"
  lazy val UseGcsIntegrationTaggedTests = Tests.Argument(TestFrameworks.ScalaTest, "-n", GcsIntegrationTestTag)
  lazy val DontUseGcsIntegrationTaggedTests = Tests.Argument(TestFrameworks.ScalaTest, "-l", GcsIntegrationTestTag)

  lazy val AwsIntegrationTestTag = "AwsTest"
  lazy val UseAwsIntegrationTaggedTests = Tests.Argument(TestFrameworks.ScalaTest, "-n", AwsIntegrationTestTag)
  lazy val DontUseAwsIntegrationTaggedTests = Tests.Argument(TestFrameworks.ScalaTest, "-l", AwsIntegrationTestTag)

  lazy val DbmsTestTag = "DbmsTest"
  lazy val UseDbmsTaggedTests = Tests.Argument(TestFrameworks.ScalaTest, "-n", DbmsTestTag)
  lazy val DontUseDbmsTaggedTests = Tests.Argument(TestFrameworks.ScalaTest, "-l", DbmsTestTag)

  lazy val TestReportArgs = Tests.Argument(TestFrameworks.ScalaTest, "-oDSI", "-h", "target/test-reports")

  /*
  The arguments that will be added to the default test config, but removed from all other configs.
  `sbt coverage test` adds other arguments added to generate the coverage reports.
  Tracking the arguments we add to the default allows one to later remove them when building up other configurations.
 */
  lazy val defaultExcludeTests = Seq(DontUseDockerTaggedTests, DontUseCromwellIntegrationTaggedTests,
    DontUseDbmsTaggedTests, DontUseGcsIntegrationTaggedTests, DontUseAwsIntegrationTaggedTests)

  val testSettings = List(
    libraryDependencies ++= testDependencies.map(_ % Test),
    // `test` (or `assembly`) - Run all tests, except docker and integration and DBMS
    testOptions in Test ++= Seq(TestReportArgs) ++ defaultExcludeTests,
    // `alltests:test` - Run all tests
    testOptions in AllTests := (testOptions in Test).value.diff(defaultExcludeTests),
    // `docker:test` - Run only the docker tests.
    testOptions in DockerTest := (testOptions in AllTests).value ++ Seq(UseDockerTaggedTests),
    // `nodocker:test` - Run all tests, except docker
    testOptions in NoDockerTest := (testOptions in AllTests).value ++ Seq(DontUseDockerTaggedTests),
    // `integration:test` - Run only integration tests
    testOptions in CromwellIntegrationTest := (testOptions in AllTests).value ++
      Seq(UseCromwellIntegrationTaggedTests, UseGcsIntegrationTaggedTests, UseAwsIntegrationTaggedTests),
    // `nointegration:test` - Run all tests, except integration
    testOptions in CromwellNoIntegrationTest := (testOptions in AllTests).value ++
      Seq(DontUseCromwellIntegrationTaggedTests, DontUseGcsIntegrationTaggedTests, DontUseAwsIntegrationTaggedTests),
    // `dbms:test` - Run database management tests.
    testOptions in DbmsTest := (testOptions in AllTests).value ++ Seq(UseDbmsTaggedTests),
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
      .configs(DockerTest).settings(inConfig(DockerTest)(Defaults.testTasks): _*)
      .configs(NoDockerTest).settings(inConfig(NoDockerTest)(Defaults.testTasks): _*)
      .configs(CromwellIntegrationTest).settings(inConfig(CromwellIntegrationTest)(Defaults.testTasks): _*)
      .configs(CromwellBenchmarkTest).settings(inConfig(CromwellBenchmarkTest)(Defaults.testTasks): _*)
      .configs(CromwellNoIntegrationTest).settings(inConfig(CromwellNoIntegrationTest)(Defaults.testTasks): _*)
      .configs(DbmsTest).settings(inConfig(DbmsTest)(Defaults.testTasks): _*)
  }

  def addIntegrationTestSettings(project: Project) = {
    project
      .settings(integrationTestSettings)
      .configs(IntegrationTest)
  }

}
