import Dependencies._
import Settings._

// Libraries

lazy val common = project
  .withLibrarySettings("cromwell-common", commonDependencies, crossCompile = true)

lazy val wom = project
  .withLibrarySettings("cromwell-wom", womDependencies, crossCompile = true)
  .dependsOn(common)

lazy val wdlRoot = Path("wdl")

lazy val wdlModelRoot = wdlRoot / "model"

lazy val wdlSharedModel = (project in wdlModelRoot / "shared")
  .withLibrarySettings("cromwell-wdl-model-core", wdlDependencies, crossCompile = true)
  .dependsOn(wom)

lazy val wdlModelDraft2 = (project in wdlModelRoot / "draft2")
  .withLibrarySettings("cromwell-wdl-model-draft2", wdlDependencies, crossCompile = true)
  .dependsOn(wdlSharedModel)

lazy val wdlTransformsRoot = wdlRoot / "transforms"

lazy val wdlSharedTransforms = (project in wdlTransformsRoot / "shared")
  .withLibrarySettings("cromwell-wdl-transforms-shared", wdlDependencies, crossCompile = true)
  .dependsOn(wdlSharedModel)
  .dependsOn(wom)

lazy val wdlTransformsDraft2 = (project in wdlTransformsRoot / "draft2")
  .withLibrarySettings("cromwell-wdl-transforms-draft2", wdlDependencies, crossCompile = true)
  .dependsOn(wdlSharedTransforms)
  .dependsOn(wdlModelDraft2)

lazy val cwl = project
  .withLibrarySettings("cromwell-cwl", cwlDependencies, crossCompile = true)
  .dependsOn(wom)
  .dependsOn(wom % "test->test")

lazy val cwlEncoder = project
  .withLibrarySettings("cromwell-cwl-encoder", crossCompile = true)
  .dependsOn(cwl)

lazy val core = project
  .withLibrarySettings("cromwell-core", coreDependencies)
  .dependsOn(wom)
  .dependsOn(wom % "test->test")

lazy val cloudSupport = project
  .withLibrarySettings("cromwell-cloud-support", cloudSupportDependencies)
  .dependsOn(common)

lazy val gcsFileSystem = (project in file("filesystems/gcs"))
  .withLibrarySettings("cromwell-gcsfilesystem")
  .dependsOn(core)
  .dependsOn(cloudSupport)
  .dependsOn(core % "test->test")
  .dependsOn(cloudSupport % "test->test")

lazy val ossFileSystem = (project in file("filesystems/oss"))
  .withLibrarySettings("cromwell-ossFileSystem", ossFileSystemDependencies)
  .dependsOn(core)
  .dependsOn(core % "test->test")

lazy val databaseSql = (project in file("database/sql"))
  .withLibrarySettings("cromwell-database-sql", databaseSqlDependencies)

lazy val databaseMigration = (project in file("database/migration"))
  .withLibrarySettings("cromwell-database-migration", databaseMigrationDependencies)
  .dependsOn(core)
  .dependsOn(wdlModelDraft2)
  .dependsOn(wdlTransformsDraft2)

lazy val dockerHashing = project
  .withLibrarySettings("cromwell-docker-hashing")
  .dependsOn(core)
  .dependsOn(core % "test->test")

lazy val cromwellApiClient = project
  .withLibrarySettings("cromwell-api-client", cromwellApiClientDependencies)

lazy val centaur = project
  .withLibrarySettings("centaur", centaurDependencies, integrationTests = true)
  .dependsOn(cloudSupport)
  .dependsOn(cromwellApiClient)

lazy val services = project
  .withLibrarySettings("cromwell-services")
  .dependsOn(databaseSql)
  .dependsOn(databaseMigration)
  .dependsOn(cloudSupport)
  .dependsOn(dockerHashing)
  .dependsOn(core % "test->test")

lazy val backendRoot = Path("supportedBackends")

lazy val backend = project
  .withLibrarySettings("cromwell-backend", backendDependencies, backendSettings)
  .dependsOn(services)
  .dependsOn(core % "test->test")

lazy val sfsBackend = (project in backendRoot / "sfs")
  .withLibrarySettings("cromwell-sfs-backend")
  .dependsOn(backend)
  .dependsOn(gcsFileSystem)
  .dependsOn(backend % "test->test")
  .dependsOn(services % "test->test")

lazy val tesBackend = (project in backendRoot / "tes")
  .withLibrarySettings("cromwell-tes-backend", tesBackendDependencies)
  .dependsOn(sfsBackend)
  .dependsOn(backend % "test->test")

lazy val sparkBackend = (project in backendRoot / "spark")
  .withLibrarySettings("cromwell-spark-backend", sparkBackendDependencies)
  .dependsOn(sfsBackend)
  .dependsOn(backend % "test->test")

lazy val jesBackend = (project in backendRoot / "jes")
  .withLibrarySettings("cromwell-jes-backend")
  .dependsOn(backend)
  .dependsOn(gcsFileSystem)
  .dependsOn(backend % "test->test")
  .dependsOn(gcsFileSystem % "test->test")
  .dependsOn(services % "test->test")

lazy val bcsBackend = (project in backendRoot / "bcs")
  .withLibrarySettings("cromwell-bcs-backend", bcsBackendDependencies)
  .dependsOn(backend)
  .dependsOn(ossFileSystem)
  .dependsOn(gcsFileSystem)
  .dependsOn(core % "test->test")
  .dependsOn(backend % "test->test")
  .dependsOn(ossFileSystem % "test->test")
  .dependsOn(services % "test->test")

lazy val engine = project
  .withLibrarySettings("cromwell-engine", engineDependencies, engineSettings)
  .dependsOn(backend)
  .dependsOn(ossFileSystem)
  .dependsOn(gcsFileSystem)
  .dependsOn(languageFactoryCore)
  .dependsOn(core % "test->test")
  .dependsOn(backend % "test->test")
  // In the future we may have a dedicated test backend like the `TestLocalAsyncJobExecutionActor`.
  // For now, all the engine tests run on the "Local" backend, an implementation of an impl.sfs.config backend.
  .dependsOn(sfsBackend % "test->compile")
  .dependsOn(gcsFileSystem % "test->test")
  .dependsOn(ossFileSystem % "test->test")

// Executables

lazy val centaurCwlRunner = project
  .withExecutableSettings("centaur-cwl-runner", centaurCwlRunnerDependencies, buildDocker = false)
  .dependsOn(cwl)
  .dependsOn(centaur)
  .dependsOn(gcsFileSystem)

lazy val womtool = project
  .withExecutableSettings("womtool", womtoolDependencies)
  .dependsOn(wdlModelDraft2)
  .dependsOn(wdlTransformsDraft2)
  .dependsOn(cwl)
  .dependsOn(wom % "test->test")

lazy val languageFactoryRoot = Path("languageFactories")

lazy val languageFactoryCore = (project in languageFactoryRoot / "language-factory-core")
  .withExecutableSettings("language-factory-core")
  .dependsOn(core)

lazy val wdlDraft2LanguageFactory = (project in languageFactoryRoot / "wdl-draft2")
  .withExecutableSettings("wdl-draft2")
  .dependsOn(languageFactoryCore)
  .dependsOn(wdlModelDraft2)
  .dependsOn(wdlTransformsDraft2)

lazy val cwlV1_0LanguageFactory = (project in languageFactoryRoot / "cwl-v1-0")
  .withExecutableSettings("cwl-v1-0")
  .dependsOn(languageFactoryCore)
  .dependsOn(cwl)

lazy val root = (project in file("."))
  .withExecutableSettings("cromwell", rootDependencies, rootSettings)
  // Next level of projects to include in the fat jar (their dependsOn will be transitively included)
  // Must 'dependsOn' any projects which might be loaded at runtime from the fat jar
  // NB: even if you 'dependsOn' a project, you still need to aggregate it (c'est la vie...)
  .dependsOn(engine)
  .dependsOn(jesBackend)
  .dependsOn(bcsBackend)
  .dependsOn(tesBackend)
  .dependsOn(sparkBackend)
  .dependsOn(cromwellApiClient)
  .dependsOn(wdlDraft2LanguageFactory)
  .dependsOn(cwlV1_0LanguageFactory)
  .dependsOn(engine % "test->test")
  // Full list of all sub-projects to build with the root (ex: include in `sbt test`)
  .aggregate(backend)
  .aggregate(centaur)
  .aggregate(centaurCwlRunner)
  .aggregate(common)
  .aggregate(core)
  .aggregate(cloudSupport)
  .aggregate(cromwellApiClient)
  .aggregate(cwl)
  .aggregate(databaseMigration)
  .aggregate(databaseSql)
  .aggregate(dockerHashing)
  .aggregate(engine)
  .aggregate(gcsFileSystem)
  .aggregate(ossFileSystem)
  .aggregate(jesBackend)
  .aggregate(services)
  .aggregate(bcsBackend)
  .aggregate(sfsBackend)
  .aggregate(sparkBackend)
  .aggregate(tesBackend)
  .dependsOn(wdlSharedModel)
  .dependsOn(wdlSharedTransforms)
  .dependsOn(wdlModelDraft2)
  .dependsOn(wdlTransformsDraft2)
  .aggregate(wom)
  .aggregate(womtool)
  .aggregate(languageFactoryCore)
  .aggregate(wdlDraft2LanguageFactory)
  .aggregate(cwlV1_0LanguageFactory)
  // TODO: See comment in plugins.sbt regarding SBT 1.x
  .enablePlugins(CrossPerProjectPlugin)
