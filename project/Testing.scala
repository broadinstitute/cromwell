import sbt.Keys._
import sbt._

object Testing {
  lazy val AllTests = config("alltests") extend Test
  lazy val NoTests = config("notests") extend Test
  lazy val DockerTest = config("docker") extend Test
  lazy val NoDockerTest = config("nodocker") extend Test
  lazy val CromwellIntegrationTest = config("integration") extend Test
  lazy val CromwellNoIntegrationTest = config("nointegration") extend Test
  lazy val DbmsTest = config("dbms") extend Test

  lazy val DockerTestTag = "DockerTest"
  lazy val UseDockerTaggedTests = Tests.Argument("-n", DockerTestTag)
  lazy val DontUseDockerTaggedTests = Tests.Argument("-l", DockerTestTag)

  lazy val CromwellIntegrationTestTag = "CromwellIntegrationTest"
  lazy val UseCromwellIntegrationTaggedTests = Tests.Argument("-n", CromwellIntegrationTestTag)
  lazy val DontUseCromwellIntegrationTaggedTests = Tests.Argument("-l", CromwellIntegrationTestTag)

  lazy val DbmsTestTag = "DbmsTest"
  lazy val UseDbmsTaggedTests = Tests.Argument("-n", DbmsTestTag)
  lazy val DontUseDbmsTaggedTests = Tests.Argument("-l", DbmsTestTag)

  /*
  The arguments that will be added to the default test config, but removed from all other configs.
  `sbt coverage test` adds other arguments added to generate the coverage reports.
  Tracking the arguments we add to the default allows one to later remove them when building up other configurations.
 */
  lazy val defaultTestArgs = Seq(DontUseDockerTaggedTests, DontUseCromwellIntegrationTaggedTests, DontUseDbmsTaggedTests)

  val testSettings = List(
    // `test` (or `assembly`) - Run all tests, except docker and integration and DBMS
    testOptions in Test ++= defaultTestArgs,
    // `alltests:test` - Run all tests
    testOptions in AllTests := (testOptions in Test).value.diff(defaultTestArgs),
    // `notests:test` - Run no tests
    testOptions in NoTests := (testOptions in Test).value ++ Seq(Tests.Filter(_ => false)),
    // `docker:test` - Run only the docker tests.
    testOptions in DockerTest := (testOptions in AllTests).value ++ Seq(UseDockerTaggedTests),
    // `nodocker:test` - Run all tests, except docker
    testOptions in NoDockerTest := (testOptions in AllTests).value ++ Seq(DontUseDockerTaggedTests),
    // `integration:test` - Run only integration tests
    testOptions in CromwellIntegrationTest := (testOptions in AllTests).value ++ Seq(UseCromwellIntegrationTaggedTests),
    // `nointegration:test` - Run all tests, except integration
    testOptions in CromwellNoIntegrationTest := (testOptions in AllTests).value ++ Seq(DontUseCromwellIntegrationTaggedTests),
    // `dbms:test` - Run database management tests.
    testOptions in DbmsTest := (testOptions in AllTests).value ++ Seq(UseDbmsTaggedTests)
  )
  /*
    TODO: This syntax of test in (NoTests, assembly) isn't correct

    Trying to get:

      sbt notests:assembly

    To be the same as:

      sbt 'set test in assembly := {}' assembly

    For now, one must use the more verbose command line until someone can crack sbt's custom configs/tasks/scopes for the
    assembly plugin:

      http://www.scala-sbt.org/0.13/tutorial/Scopes.html
   */
  //test in (NoTests, assembly) := {}

  // Also tried
  //test in (NoTests, assemblyPackageDependency) := {}
  //test in (NoTests, assemblyPackageScala) := {}
}
