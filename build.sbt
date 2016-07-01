import Settings._
import Testing._

lazy val core = (project in file("core"))
  .settings(coreSettings:_*)
  .withTestSettings

lazy val gcsfilesystem = (project in file("filesystems/gcs"))
  .settings(gcsFileSystemSettings:_*)
  .withTestSettings

lazy val database = (project in file("database"))
  .settings(databaseSettings:_*)
  .withTestSettings

lazy val services = (project in file("services"))
  .settings(servicesSettings:_*)
  .dependsOn(core % "test->test;compile->compile")
  .dependsOn(database % "test->test;compile->compile") // Assuming services impl in services
  .withTestSettings

lazy val backendRoot = Path("supportedBackends")

lazy val backend = (project in file("backend"))
  .dependsOn(core % "test->test;compile->compile")
  .dependsOn(services % "test->test;compile->compile")
  .settings(backendSettings:_*)
  .withTestSettings

lazy val localBackend = (project in backendRoot / "local")
  .dependsOn(backend % "test->test;compile->compile")
  .settings(localBackendSettings:_*)
  .withTestSettings

lazy val htCondorBackend = (project in backendRoot / "htcondor")
  .settings(htCondorBackendSettings:_*)
  .dependsOn(backend % "test->test;compile->compile")
  .withTestSettings

lazy val sparkBackend = (project in backendRoot / "spark")
  .settings(sparkBackendSettings:_*)
  .dependsOn(backend % "test->test;compile->compile")
  .withTestSettings

lazy val sgeBackend = (project in backendRoot / "sge")
  .settings(sgeBackendSettings:_*)
  .dependsOn(backend % "test->test;compile->compile")
  .withTestSettings

lazy val jesBackend = (project in backendRoot / "jes")
  .settings(jesBackendSettings:_*)
  .dependsOn(backend % "test->test;compile->compile")
  .dependsOn(gcsfilesystem % "test->test;compile->compile")
  .withTestSettings

//TODO: remove jesBackend once refactoring has finished.
lazy val engine = (project in file("engine"))
  .settings(engineSettings: _*)
  .dependsOn(core % "test->test;compile->compile")
  .dependsOn(services % "test->test;compile->compile")
  .dependsOn(backend % "test->test;compile->compile")
  .dependsOn(jesBackend % "test->test;compile->compile")
  .dependsOn(localBackend % "test->test;compile->compile")
  .dependsOn(htCondorBackend % "test->test;compile->compile")
  .dependsOn(database % "test->test;compile->compile")
  .withTestSettings

lazy val root = (project in file("."))
  .settings(rootSettings: _*)
  .dependsOn(core % "test->test;compile->compile")
  .dependsOn(engine % "test->test;compile->compile")
  .aggregate(core, database, backend, engine, localBackend, sgeBackend, jesBackend, htCondorBackend, sparkBackend, gcsfilesystem)
  .withTestSettings
