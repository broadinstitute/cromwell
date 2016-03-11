import Settings._
import Testing._

lazy val core = (project in file("core")).settings(coreSettings:_*)

lazy val backend = (project in file("backend"))
  .dependsOn(core % "test->test;compile->compile")
  .settings(backendSettings:_*)

lazy val engine = (project in file("engine"))
  .dependsOn(core % "test->test;compile->compile", backend)
  .settings(engineSettings:_*)
  .configs(AllTests).settings(inConfig(AllTests)(Defaults.testTasks): _*)
  .configs(NoTests).settings(inConfig(NoTests)(Defaults.testTasks): _*)
  .configs(DockerTest).settings(inConfig(DockerTest)(Defaults.testTasks): _*)
  .configs(NoDockerTest).settings(inConfig(NoDockerTest)(Defaults.testTasks): _*)
  .configs(CromwellIntegrationTest).settings(inConfig(CromwellIntegrationTest)(Defaults.testTasks): _*)
  .configs(CromwellNoIntegrationTest).settings(inConfig(CromwellNoIntegrationTest)(Defaults.testTasks): _*)
  .configs(DbmsTest).settings(inConfig(DbmsTest)(Defaults.testTasks): _*)

lazy val root = (project in file("."))
  .dependsOn(engine % "test->test;compile->compile", core % "test->test;compile->compile", backend)
  .aggregate(core, backend, engine)
  .settings(rootSettings:_*)
  .configs(AllTests).settings(inConfig(AllTests)(Defaults.testTasks): _*)
  .configs(NoTests).settings(inConfig(NoTests)(Defaults.testTasks): _*)
  .configs(DockerTest).settings(inConfig(DockerTest)(Defaults.testTasks): _*)
  .configs(NoDockerTest).settings(inConfig(NoDockerTest)(Defaults.testTasks): _*)
  .configs(CromwellIntegrationTest).settings(inConfig(CromwellIntegrationTest)(Defaults.testTasks): _*)
  .configs(CromwellNoIntegrationTest).settings(inConfig(CromwellNoIntegrationTest)(Defaults.testTasks): _*)
  .configs(DbmsTest).settings(inConfig(DbmsTest)(Defaults.testTasks): _*)


