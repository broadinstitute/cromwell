import Dependencies._
import Settings._

// Libraries

lazy val common = project
  .withLibrarySettings("cromwell-common", commonDependencies)

lazy val wom = project
  .withLibrarySettings("cromwell-wom", womDependencies)
  .dependsOn(common)

lazy val wdlRoot = Path("wdl")

lazy val wdlModelRoot = wdlRoot / "model"

lazy val wdlSharedModel = (project in wdlModelRoot / "shared")
  .withLibrarySettings("cromwell-wdl-model-core", wdlDependencies)
  .dependsOn(wom)

lazy val wdlModelDraft2 = (project in wdlModelRoot / "draft2")
  .withLibrarySettings("cromwell-wdl-model-draft2", wdlDependencies)
  .dependsOn(wdlSharedTransforms)
  .dependsOn(wdlSharedModel)
  .dependsOn(wom % "test->test")

lazy val wdlModelDraft3 = (project in wdlModelRoot / "draft3")
  .withLibrarySettings("cromwell-wdl-model-draft3")
  .dependsOn(wdlSharedModel)

lazy val wdlTransformsRoot = wdlRoot / "transforms"

lazy val wdlSharedTransforms = (project in wdlTransformsRoot / "shared")
  .withLibrarySettings("cromwell-wdl-transforms-shared", wdlDependencies)
  .dependsOn(wdlSharedModel)
  .dependsOn(wom)

lazy val wdlTransformsDraft2 = (project in wdlTransformsRoot / "draft2")
  .withLibrarySettings("cromwell-wdl-transforms-draft2", wdlDependencies)
  .dependsOn(wdlSharedTransforms)
  .dependsOn(wdlModelDraft2)

lazy val wdlTransformsDraft3 = (project in wdlTransformsRoot / "draft3")
  .withLibrarySettings("cromwell-wdl-transforms-draft3", wdlDependencies)
  .dependsOn(wdlSharedTransforms)
  .dependsOn(wdlModelDraft3)
  .dependsOn(languageFactoryCore)
  .dependsOn(common % "test->test")
  .dependsOn(wom % "test->test")

lazy val cwl = project
  .withLibrarySettings("cromwell-cwl", cwlDependencies)
  .dependsOn(wom)
  .dependsOn(wom % "test->test")

lazy val cwlEncoder = project
  .withLibrarySettings("cromwell-cwl-encoder")
  .dependsOn(cwl)

lazy val core = project
  .withLibrarySettings("cromwell-core", coreDependencies)
  .dependsOn(wom)
  .dependsOn(wom % "test->test")

lazy val cloudSupport = project
  .withLibrarySettings("cromwell-cloud-support", cloudSupportDependencies)
  .dependsOn(common)

lazy val awsS3FileSystem = (project in file("filesystems/s3"))
  .withLibrarySettings("cromwell-aws-s3filesystem")
  .dependsOn(core)
  .dependsOn(cloudSupport)
  .dependsOn(core % "test->test")
  .dependsOn(cloudSupport % "test->test")

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

lazy val googlePipelinesCommon = (project in backendRoot / "google" / "pipelines" / "common")
  .withLibrarySettings("cromwell-pipelines-common")
  .dependsOn(backend)
  .dependsOn(gcsFileSystem)
  .dependsOn(backend % "test->test")
  .dependsOn(gcsFileSystem % "test->test")
  .dependsOn(services % "test->test")

lazy val googlePipelinesV1Alpha2 = (project in backendRoot / "google" / "pipelines" / "v1alpha2")
  .withLibrarySettings("cromwell-pipelines-v1-backend")
  .dependsOn(googlePipelinesCommon)
  .dependsOn(googlePipelinesCommon % "test->test")
  .dependsOn(core % "test->test")

lazy val googlePipelinesV2Alpha1 = (project in backendRoot / "google" / "pipelines" / "v2alpha1")
  .withLibrarySettings("cromwell-pipelines-v2-backend")
  .dependsOn(googlePipelinesCommon)
  .dependsOn(googlePipelinesCommon % "test->test")
  .dependsOn(core % "test->test")

// Legacy, inherits all its code from googlePipelinesV1Alpha2
lazy val jesBackend = (project in backendRoot / "jes")
  .withLibrarySettings("cromwell-jes-backend")
  .dependsOn(googlePipelinesV1Alpha2)

lazy val awsBackend = (project in backendRoot / "aws")
  .withLibrarySettings("cromwell-aws-backend")
  .dependsOn(backend)
  .dependsOn(awsS3FileSystem)
  .dependsOn(backend % "test->test")
  .dependsOn(awsS3FileSystem % "test->test")
  .dependsOn(services % "test->test")

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
  .dependsOn(wdlDraft2LanguageFactory)
  .dependsOn(wdlDraft3LanguageFactory)
  .dependsOn(cwlV1_0LanguageFactory)
  .dependsOn(wom % "test->test")

lazy val cromiam = (project in file("CromIAM")) // TODO: git mv CromIAM to a canonical lowercased name
  .withExecutableSettings("cromiam", cromiamDependencies, cromiamSettings)
  .dependsOn(common)
  .dependsOn(cromwellApiClient)

lazy val languageFactoryRoot = Path("languageFactories")

lazy val languageFactoryCore = (project in languageFactoryRoot / "language-factory-core")
  .withLibrarySettings("language-factory-core", languageFactoryDependencies)
  .dependsOn(core)

lazy val wdlDraft2LanguageFactory = (project in languageFactoryRoot / "wdl-draft2")
  .withLibrarySettings("wdl-draft2", draft2LanguageFactoryDependencies)
  .dependsOn(languageFactoryCore)
  .dependsOn(wdlModelDraft2)
  .dependsOn(wdlTransformsDraft2)

lazy val wdlDraft3LanguageFactory = (project in languageFactoryRoot / "wdl-draft3")
  .withLibrarySettings("wdl-draft3")
  .dependsOn(languageFactoryCore)
  .dependsOn(wdlModelDraft3)
  .dependsOn(wdlTransformsDraft3)

lazy val cwlV1_0LanguageFactory = (project in languageFactoryRoot / "cwl-v1-0")
  .withLibrarySettings("cwl-v1-0")
  .dependsOn(languageFactoryCore)
  .dependsOn(cwl)

lazy val server = project
  .withExecutableSettings("cromwell", serverDependencies)
  .dependsOn(engine)
  .dependsOn(googlePipelinesV1Alpha2)
  .dependsOn(googlePipelinesV2Alpha1)
  .dependsOn(jesBackend)
  .dependsOn(bcsBackend)
  .dependsOn(awsBackend)
  .dependsOn(tesBackend)
  .dependsOn(sparkBackend)
  .dependsOn(cromwellApiClient)
  .dependsOn(wdlDraft2LanguageFactory)
  .dependsOn(wdlDraft3LanguageFactory)
  .dependsOn(cwlV1_0LanguageFactory)
  .dependsOn(engine % "test->test")

lazy val root = (project in file("."))
  .withRootSettings()
  // Full list of all sub-projects to build with the root (ex: include in `sbt test`)
  .aggregate(backend)
  .aggregate(bcsBackend)
  .aggregate(centaur)
  .aggregate(centaurCwlRunner)
  .aggregate(cloudSupport)
  .aggregate(common)
  .aggregate(core)
  .aggregate(cromiam)
  .aggregate(cromwellApiClient)
  .aggregate(cwl)
  .aggregate(cwlV1_0LanguageFactory)
  .aggregate(databaseMigration)
  .aggregate(databaseSql)
  .aggregate(dockerHashing)
  .aggregate(engine)
  .aggregate(awsS3FileSystem)
  .aggregate(gcsFileSystem)
  .aggregate(googlePipelinesCommon)
  .aggregate(googlePipelinesV1Alpha2)
  .aggregate(googlePipelinesV2Alpha1)
  .aggregate(jesBackend)
  .aggregate(languageFactoryCore)
  .aggregate(awsBackend)
  .aggregate(ossFileSystem)
  .aggregate(server)
  .aggregate(services)
  .aggregate(sfsBackend)
  .aggregate(sparkBackend)
  .aggregate(tesBackend)
  .aggregate(wdlDraft2LanguageFactory)
  .aggregate(wdlDraft3LanguageFactory)
  .aggregate(wdlModelDraft2)
  .aggregate(wdlModelDraft3)
  .aggregate(wdlSharedModel)
  .aggregate(wdlSharedTransforms)
  .aggregate(wdlTransformsDraft2)
  .aggregate(wdlTransformsDraft3)
  .aggregate(wom)
  .aggregate(womtool)
