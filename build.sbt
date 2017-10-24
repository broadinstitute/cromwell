import Settings._
import Testing._

lazy val common = (project in file("common"))
  .settings(commonSettings:_*)
  .withTestSettings

lazy val wom = (project in file("wom"))
  .settings(womSettings:_*)
  .dependsOn(common)
  .withTestSettings

lazy val wdl = (project in file("wdl"))
  .settings(wdlSettings:_*)
  .dependsOn(wom)
  .withTestSettings

lazy val cwl = (project in file("cwl"))
  .settings(cwlSettings:_*)
  .dependsOn(wom)
  .dependsOn(wom % "test->test")
  .withTestSettings

lazy val womtool = (project in file("womtool"))
  .settings(womtoolSettings:_ *)
  .dependsOn(wdl)
  .dependsOn(cwl)
  .dependsOn(wom % "test->test")
  .withTestSettings

lazy val core = (project in file("core"))
  .settings(coreSettings:_*)
  .dependsOn(wom)
  .dependsOn(wom % "test->test")
  .withTestSettings

lazy val gcsFileSystem = (project in file("filesystems/gcs"))
  .settings(gcsFileSystemSettings:_*)
  .withTestSettings
  .dependsOn(core)
  .dependsOn(core % "test->test")
  .dependsOn(cloudSupport)
  .dependsOn(cloudSupport % "test->test")

lazy val databaseSql = (project in file("database/sql"))
  .settings(databaseSqlSettings:_*)
  .withTestSettings

lazy val databaseMigration = (project in file("database/migration"))
  .settings(databaseMigrationSettings: _*)
  .dependsOn(core)
  .dependsOn(wdl)
  .withTestSettings

lazy val dockerHashing = (project in file("dockerHashing"))
  .settings(dockerHashingSettings: _*)
  .dependsOn(core)
  .dependsOn(core % "test->test")
  .withTestSettings

lazy val cromwellApiClient = (project in file("cromwellApiClient"))
  .settings(cromwellApiClientSettings: _*)
  .withTestSettings

lazy val services = (project in file("services"))
  .settings(servicesSettings:_*)
  .withTestSettings
  .dependsOn(core)
  .dependsOn(databaseSql)
  .dependsOn(databaseMigration)
  .dependsOn(cloudSupport)
  .dependsOn(dockerHashing)
  .dependsOn(core % "test->test")

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
  .dependsOn(services % "test->test")

lazy val tesBackend = (project in backendRoot / "tes")
  .settings(tesBackendSettings:_*)
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
  .dependsOn(cloudSupport)
  .dependsOn(backend % "test->test")
  .dependsOn(gcsFileSystem % "test->test")
  .dependsOn(services % "test->test")

lazy val engine = (project in file("engine"))
  .settings(engineSettings: _*)
  .withTestSettings
  .dependsOn(core)
  .dependsOn(dockerHashing)
  .dependsOn(services)
  .dependsOn(backend)
  .dependsOn(gcsFileSystem)
  .dependsOn(wdl)
  .dependsOn(cwl)
  .dependsOn(core % "test->test")
  .dependsOn(backend % "test->test")
  // In the future we may have a dedicated test backend like the `TestLocalAsyncJobExecutionActor`.
  // For now, all the engine tests run on the "Local" backend, an implementation of an impl.sfs.config backend.
  .dependsOn(sfsBackend % "test->compile")
  .dependsOn(gcsFileSystem % "test->test")

lazy val cloudSupport = (project in file("cloudSupport"))
    .settings(cloudSupportSettings: _*)
    .withTestSettings
    .dependsOn(core)
    .dependsOn(core % "test->test")

lazy val root = (project in file("."))
  .settings(rootSettings: _*)
  .enablePlugins(DockerPlugin)
  .withTestSettings
  // Full list of all sub-projects to build with the root (ex: include in `sbt test`)
  .aggregate(common)
  .aggregate(wom)
  .aggregate(wdl)
  .aggregate(cwl)
  .aggregate(womtool)
  .aggregate(core)
  .aggregate(dockerHashing)
  .aggregate(gcsFileSystem)
  .aggregate(databaseSql)
  .aggregate(databaseMigration)
  .aggregate(services)
  .aggregate(backend)
  .aggregate(sfsBackend)
  .aggregate(sparkBackend)
  .aggregate(jesBackend)
  .aggregate(tesBackend)
  .aggregate(engine)
  .aggregate(cromwellApiClient)
  // Next level of projects to include in the fat jar (their dependsOn will be transitively included)
  .dependsOn(engine)
  .dependsOn(jesBackend)
  .dependsOn(tesBackend)
  .dependsOn(sparkBackend)
  .dependsOn(wdl)
  .dependsOn(cwl)
  .dependsOn(womtool)
  // Dependencies for tests
  .dependsOn(engine % "test->test")

logLevel in ThisBuild := Level.Warn
