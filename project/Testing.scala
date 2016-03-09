import sbt.Keys._
import sbt._

object Testing {
  lazy val AllTests = config("alltests") extend Test
  lazy val NoTests = config("notests") extend Test
  lazy val DockerTest = config("docker") extend Test
  lazy val NoDockerTest = config("nodocker") extend Test
  lazy val CromwellIntegrationTest = config("integration") extend Test
  lazy val CromwellNoIntegrationTest = config("nointegration") extend Test

  /*
  The arguments that will be added to the default test config, but removed from all other configs.
  `sbt coverage test` adds other arguments added to generate the coverage reports.
  Tracking the arguments we add to the default allows one to later remove them when building up other configurations.
 */
  lazy val defaultTestArgs = Seq(Tests.Argument("-l", "DockerTest"), Tests.Argument("-l", "CromwellIntegrationTest"))

  val testSettings = List(
    // `test` (or `assembly`) - Run all tests, except docker and integration
    testOptions in Test ++= defaultTestArgs,
    // `alltests:test` - Run all tests
    testOptions in AllTests := (testOptions in Test).value.diff(defaultTestArgs),
    // `docker:test` - Run docker tests, except integration
    testOptions in DockerTest := (testOptions in Test).value.diff(defaultTestArgs) ++
      Seq(Tests.Argument("-n", "DockerTest"), Tests.Argument("-l", "CromwellIntegrationTest")),
    // `nodocker:test` - Run all tests, except docker
    testOptions in NoDockerTest := (testOptions in Test).value.diff(defaultTestArgs) ++
      Seq(Tests.Argument("-l", "DockerTest")),
    // `integration:test` - Run integration tests, except docker
    testOptions in CromwellIntegrationTest := (testOptions in Test).value.diff(defaultTestArgs) ++
      Seq(Tests.Argument("-l", "DockerTest"), Tests.Argument("-n", "CromwellIntegrationTest")),
    // `nointegration:test` - Run all tests, except integration
    testOptions in CromwellNoIntegrationTest := (testOptions in Test).value.diff(defaultTestArgs) ++
      Seq(Tests.Argument("-l", "CromwellIntegrationTest"))
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
