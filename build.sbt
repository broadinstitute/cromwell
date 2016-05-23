import Settings._
import Testing._

lazy val core = (project in file("core")).settings(coreSettings:_*)

lazy val gcsfilesystem = (project in file("filesystems/gcs"))
  .settings(gcsFileSystemSettings:_*)
  .configs(AllTests).settings(inConfig(AllTests)(Defaults.testTasks): _*)
  .configs(CromwellIntegrationTest).settings(inConfig(CromwellIntegrationTest)(Defaults.testTasks): _*)
  .configs(CromwellNoIntegrationTest).settings(inConfig(CromwellNoIntegrationTest)(Defaults.testTasks): _*)

lazy val database = (project in file("database"))
  .settings(databaseSettings:_*)
  .configs(AllTests).settings(inConfig(AllTests)(Defaults.testTasks): _*)
  .configs(CromwellIntegrationTest).settings(inConfig(CromwellIntegrationTest)(Defaults.testTasks): _*)
  .configs(CromwellNoIntegrationTest).settings(inConfig(CromwellNoIntegrationTest)(Defaults.testTasks): _*)
  .configs(DbmsTest).settings(inConfig(DbmsTest)(Defaults.testTasks): _*)

lazy val services = (project in file("services"))
  .dependsOn(core % "test->test;compile->compile")
  .dependsOn(database % "test->test;compile->compile") // Assuming services impl in services
  .settings(servicesSettings:_*)
  .configs(AllTests).settings(inConfig(AllTests)(Defaults.testTasks): _*)
  .configs(CromwellIntegrationTest).settings(inConfig(CromwellIntegrationTest)(Defaults.testTasks): _*)
  .configs(CromwellNoIntegrationTest).settings(inConfig(CromwellNoIntegrationTest)(Defaults.testTasks): _*)

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
  .dependsOn(backend % "test->test;compile->compile")
  .dependsOn(services % "test->test;compile->compile")
  .dependsOn(gcsfilesystem % "test->test;compile->compile")
  .settings(jesBackendSettings:_*)

//TODO: remove jesBackend once refactoring has finished.
lazy val engine = (project in file("engine"))
  .dependsOn(core % "test->test;compile->compile",
    services % "test->test;compile->compile",
    backend % "test->test;compile->compile",
    jesBackend % "test->test;compile->compile",
    localBackend % "test->test;compile->compile",
    database % "test->test;compile->compile")
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


