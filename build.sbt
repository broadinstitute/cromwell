import Settings._
import Testing._

lazy val core = (project in file("core")).settings(coreSettings:_*)

lazy val gcsfilesystem = (project in file("filesystems/gcs"))
  .settings(gcsFileSystemSettings:_*)

lazy val backendRoot = Path("supportedBackends")

lazy val backend = (project in file("backend"))
  .dependsOn(core % "test->test;compile->compile")
  .settings(backendSettings:_*)

lazy val localBackend = (project in backendRoot / "local")
  .dependsOn(backend % "test->test;compile->compile")
  .settings(localBackendSettings:_*)

lazy val htCondorBackend = (project in backendRoot / "htcondor")
  .dependsOn(backend % "test->test;compile->compile")
  .settings(htCondorBackendSettings:_*)

lazy val sgeBackend = (project in backendRoot / "sge")
  .dependsOn(backend % "test->test;compile->compile")
  .settings(sgeBackendSettings:_*)

lazy val jesBackend = (project in backendRoot / "jes")
  .dependsOn(backend % "test->test;compile->compile", gcsfilesystem)
  .settings(jesBackendSettings:_*)

lazy val engine = (project in file("engine"))
  .dependsOn(core % "test->test;compile->compile", gcsfilesystem % "test->test;compile->compile", backend)
  .settings(engineSettings:_*)
  .configs(AllTests).settings(inConfig(AllTests)(Defaults.testTasks): _*)
  .configs(DockerTest).settings(inConfig(DockerTest)(Defaults.testTasks): _*)
  .configs(NoDockerTest).settings(inConfig(NoDockerTest)(Defaults.testTasks): _*)
  .configs(CromwellIntegrationTest).settings(inConfig(CromwellIntegrationTest)(Defaults.testTasks): _*)
  .configs(CromwellNoIntegrationTest).settings(inConfig(CromwellNoIntegrationTest)(Defaults.testTasks): _*)
  .configs(DbmsTest).settings(inConfig(DbmsTest)(Defaults.testTasks): _*)

lazy val root = (project in file("."))
  .dependsOn(engine % "test->test;compile->compile", core % "test->test;compile->compile", backend, localBackend, sgeBackend, jesBackend, htCondorBackend, gcsfilesystem)
  .aggregate(core, backend, engine, localBackend, sgeBackend, jesBackend, htCondorBackend, gcsfilesystem)
  .settings(rootSettings:_*)
  .configs(AllTests).settings(inConfig(AllTests)(Defaults.testTasks): _*)
  .configs(DockerTest).settings(inConfig(DockerTest)(Defaults.testTasks): _*)
  .configs(NoDockerTest).settings(inConfig(NoDockerTest)(Defaults.testTasks): _*)
  .configs(CromwellIntegrationTest).settings(inConfig(CromwellIntegrationTest)(Defaults.testTasks): _*)
  .configs(CromwellNoIntegrationTest).settings(inConfig(CromwellNoIntegrationTest)(Defaults.testTasks): _*)
  .configs(DbmsTest).settings(inConfig(DbmsTest)(Defaults.testTasks): _*)


