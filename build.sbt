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
  .withTestSettings
  .dependsOn(core)
  .dependsOn(database) // Assuming services impl in services
  .dependsOn(core % "test->test")
  .dependsOn(database % "test->test")

lazy val backendRoot = Path("supportedBackends")

lazy val backend = (project in file("backend"))
  .settings(backendSettings:_*)
  .withTestSettings
  .dependsOn(core)
  .dependsOn(services)
  .dependsOn(core % "test->test")

lazy val localBackend = (project in backendRoot / "local")
  .settings(localBackendSettings:_*)
  .withTestSettings
  .dependsOn(backend)
  .dependsOn(backend % "test->test")

lazy val htCondorBackend = (project in backendRoot / "htcondor")
  .settings(htCondorBackendSettings:_*)
  .withTestSettings
  .dependsOn(backend)
  .dependsOn(backend % "test->test")

lazy val sgeBackend = (project in backendRoot / "sge")
  .settings(sgeBackendSettings:_*)
  .withTestSettings
  .dependsOn(backend)
  .dependsOn(backend % "test->test")

lazy val jesBackend = (project in backendRoot / "jes")
  .settings(jesBackendSettings:_*)
  .withTestSettings
  .dependsOn(backend)
  .dependsOn(gcsfilesystem)
  .dependsOn(backend % "test->test")
  .dependsOn(gcsfilesystem % "test->test")

//TODO: remove jesBackend once refactoring has finished.
lazy val engine = (project in file("engine"))
  .settings(engineSettings: _*)
  .withTestSettings
  .dependsOn(core)
  .dependsOn(services)
  .dependsOn(backend)
  .dependsOn(database)
  .dependsOn(gcsfilesystem)
  .dependsOn(core % "test->test")
  .dependsOn(backend % "test->test")
  .dependsOn(localBackend % "test->compile")
  .dependsOn(jesBackend % "test->compile")
  .dependsOn(gcsfilesystem % "test->test")

lazy val root = (project in file("."))
  .settings(rootSettings: _*)
  .withTestSettings
  .dependsOn(engine)
  .dependsOn(jesBackend)
  .dependsOn(localBackend)
  .dependsOn(htCondorBackend)
  .dependsOn(engine % "test->test")
  .aggregate(core, database, services, backend, engine, localBackend, sgeBackend, jesBackend, htCondorBackend, gcsfilesystem)
