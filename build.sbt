import Settings._
import Testing._

lazy val core = (project in file("core"))
  .settings(coreSettings:_*)
  .withTestSettings

lazy val gcsFileSystem = (project in file("filesystems/gcs"))
  .settings(gcsFileSystemSettings:_*)
  .withTestSettings

lazy val database = (project in file("database"))
  .settings(databaseSettings:_*)
  .dependsOn(core) // used for migration
  .dependsOn(core % "test->test") // TODO: PBE: Remove once PostMVP tag removed from test
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

lazy val sfsBackend = (project in backendRoot / "sfs")
  .settings(sfsBackendSettings:_*)
  .withTestSettings
  .dependsOn(backend)
  .dependsOn(gcsFileSystem)
  .dependsOn(backend % "test->test")

lazy val htCondorBackend = (project in backendRoot / "htcondor")
  .settings(htCondorBackendSettings:_*)
  .withTestSettings
  .dependsOn(sfsBackend)
  .dependsOn(backend % "test->test")

lazy val sparkBackend = (project in backendRoot / "spark")
  .settings(sparkBackendSettings:_*)
  .withTestSettings
  .dependsOn(sfsBackend)
  .dependsOn(backend % "test->test")

lazy val jesBackend = (project in backendRoot / "jes")
  .settings(jesBackendSettings:_*)
  .withTestSettings
  .dependsOn(backend)
  .dependsOn(gcsFileSystem)
  .dependsOn(backend % "test->test")
  .dependsOn(gcsFileSystem % "test->test")

lazy val engine = (project in file("engine"))
  .settings(engineSettings: _*)
  .withTestSettings
  .dependsOn(core)
  .dependsOn(services)
  .dependsOn(backend)
  .dependsOn(database)
  .dependsOn(gcsFileSystem)
  .dependsOn(core % "test->test")
  .dependsOn(backend % "test->test")
  // In the future we may have a dedicated test backend like the `TestLocalAsyncJobExecutionActor`.
  // For now, all the engine tests run on the "Local" backend, an implementation of an impl.sfs.config backend.
  .dependsOn(sfsBackend % "test->compile")
  .dependsOn(gcsFileSystem % "test->test")

lazy val root = (project in file("."))
  .settings(rootSettings: _*)
  .withTestSettings
  // Full list of all sub-projects to build with the root (ex: include in `sbt test`)
  .aggregate(core)
  .aggregate(gcsFileSystem)
  .aggregate(database)
  .aggregate(services)
  .aggregate(backend)
  .aggregate(sfsBackend)
  .aggregate(htCondorBackend)
  .aggregate(sparkBackend)
  .aggregate(jesBackend)
  .aggregate(engine)
  // Next level of projects to include in the fat jar (their dependsOn will be transitively included)
  .dependsOn(engine)
  .dependsOn(jesBackend)
  .dependsOn(htCondorBackend)
  .dependsOn(sparkBackend)
  // Dependencies for tests
  .dependsOn(engine % "test->test")
