import sbt._

object Dependencies {
  private val akkaHttpCirceIntegrationV = "1.39.2"
  private val akkaHttpV = "10.1.15" // (CROM-6619)
  private val akkaV = "2.5.32" // scala-steward:off (CROM-6637)
  private val ammoniteOpsV = "2.4.1"
  private val apacheHttpClientV = "4.5.13"
  private val apacheHttpClient5V = "5.3.1"
  private val awsSdkV = "2.17.265"
  // We would like to use the BOM to manage Azure SDK versions, but SBT doesn't support it.
  // https://github.com/Azure/azure-sdk-for-java/tree/main/sdk/boms/azure-sdk-bom
  // https://github.com/sbt/sbt/issues/4531
  private val azureIdentitySdkV = "1.9.1"
  private val azureIdentityExtensionsV = "1.1.4"
  private val azureCoreManagementV = "1.7.1"
  // We are using the older AppInsights 2 because we want to use the
  // logback appender to send logs. AppInsights 3 does not have a standalone
  // appender, and its auto-hoovering of logs didn't meet our needs.
  // (Specifically, the side-by-side root logger and workflow logger resulted in
  // duplicate messages in AI. See WX-1122.)
  private val azureAppInsightsLogbackV = "2.6.4"
  private val betterFilesV = "3.9.1"
  private val jsonSmartV = "2.4.10"
  private val bardClientV = "1.0.4"
  /*
  cats-effect, fs2, http4s, and sttp (also to v3) should all be upgraded at the same time to use cats-effect 3.x.
   */
  private val catsEffectV = "2.5.3" // scala-steward:off (CROM-6564)
  private val catsV = "2.7.0"
  private val circeConfigV = "0.8.0"
  private val circeGenericExtrasV = "0.14.1"
  private val circeOpticsV = "0.14.1"
  private val circeV = "0.14.1"
  private val circeYamlV = "0.14.1"
  private val commonsCodecV = "1.15" // via: https://commons.apache.org/proper/commons-codec/
  private val commonsCsvV = "1.9.0"
  private val commonsIoV = "2.11.0" // via: https://commons.apache.org/proper/commons-io/
  private val commonsLang3V = "3.12.0"
  private val commonsMathV = "3.6.1"
  private val commonNetV = "3.8.0" // via: https://commons.apache.org/proper/commons-net/
  private val commonsTextV = "1.10.0"
  private val configsV = "0.6.1"
  private val delightRhinoSandboxV = "0.0.15"
  private val diffsonSprayJsonV = "4.1.1"
  private val ficusV = "1.5.2"
  private val fs2V = "2.5.9" // scala-steward:off (CROM-6564)
  private val googleApiClientV = "2.1.4"
  private val googleCloudBigQueryV = "2.25.0"
  // latest date via: https://github.com/googleapis/google-api-java-client-services/blob/main/clients/google-api-services-cloudkms/v1.metadata.json
  private val googleCloudKmsV = "v1-rev20230421-2.0.0"
  private val googleCloudMonitoringV = "3.2.5"
  private val googleCloudNioV = "0.124.8"
  private val googleCloudStorageV = "2.17.2"
  private val googleGaxGrpcV = "2.25.0"
  // latest date via: https://mvnrepository.com/artifact/com.google.apis/google-api-services-genomics
  private val googleGenomicsServicesV2Alpha1ApiV = "v2alpha1-rev20210811-1.32.1"
  private val googleHttpClientApacheV = "2.1.2"
  private val googleHttpClientV = "1.42.3"
  private val googleCloudBatchV1 = "0.18.0"
  // latest date via: https://mvnrepository.com/artifact/com.google.apis/google-api-services-lifesciences
  private val googleLifeSciencesServicesV2BetaApiV = "v2beta-rev20220916-2.0.0"
  private val googleOauth2V = "1.5.3"
  private val googleOauthClientV = "1.33.1"
  private val googleCloudResourceManagerV = "1.17.0"
  private val grpcV = "1.54.1"
  private val guavaV = "31.0.1-jre"
  private val heterodonV = "1.0.0-beta3"
  private val hsqldbV = "2.6.1"
  private val http4sV = "0.21.31" // this release is EOL. We need to upgrade further for cats3. https://http4s.org/versions/
  private val jacksonV = "2.13.3"
  private val jakartaActivationV = "1.2.1"
  private val jakartaAnnotationV = "1.3.5"
  private val jakartaInjectV = "2.6.1"
  private val jakartaXmlBindApiV = "2.3.2"
  private val janinoV = "3.1.7"
  private val jerseyV = "2.32" // Use a jersey compatible with WSM. See notes in wsmDependencies below.
  private val jsr305V = "3.0.2"
  private val junitV = "4.13.2"
  private val kindProjectorV = "0.13.2"
  private val kittensV = "2.3.2"
  private val liquibaseV = "4.8.0"
  private val logbackV = "1.2.13"
  private val lz4JavaV = "1.8.0"
  private val mariadbV = "2.7.4"
  /*
  The StatsD reporter for DropWizard's (Code Hale's) Metrics 3.x still works with Metrics 4.x.
  Still would be great to move to Prometheus / OpenCensus
   */
  private val metrics4ScalaV = "4.2.8"
  private val metrics3StatsdV = "4.2.0"
  private val mockFtpServerV = "3.0.0"
  private val mockitoV = "3.12.4"
  private val mockserverNettyV = "5.14.0"
  private val mouseV = "1.0.11"

  /*
  Newer version 8.0.29 fails `Control characters should work with metadata` Centaur tests, has charset changes mentioned in release notes
  https://dev.mysql.com/doc/relnotes/mysql/8.0/en/news-8-0-29.html#mysqld-8-0-29-charset
   */
  private val mysqlV = "8.0.28"
  private val nettyV = "4.1.72.Final"
  private val owlApiV = "5.1.19"
  private val pact4sV = "0.9.0"
  private val postgresV = "42.4.4"
  private val pprintV = "0.7.3"
  private val rdf4jV = "3.7.1"
  private val re2jV = "1.6"
  private val refinedV = "0.10.1"
  private val rhinoV = "1.7.14"

  private val scalaCollectionCompatV = "2.5.0"
  private val scalaGraphV = "1.13.1"
  private val scalaLoggingV = "3.9.5"
  private val scalaPoolV = "0.4.3"
  private val scalacticV = "3.2.13"
  private val scalameterV = "0.21"
  private val scalatestV = "3.2.15"
  private val scalatestScalacheckV = scalatestV + ".0"
  private val scoptV = "4.1.0"
  private val sentryLogbackV = "5.7.4"
  private val shapelessV = "2.3.9"
  private val simulacrumV = "1.0.1"
  private val slf4jV = "1.7.32"
  private val slickCatsV = "0.10.4"
  /* If you're about to update our Slick version:
    * Consider checking whether the new Slick version passes tests with upserts enabled (eg KeyValueDatabaseSpec)
    *
    * Current version 3.3.2-2076hotfix was built locally from https://github.com/grsterin/slick/tree/v3.3.2-2076hotfix
    * and manually uploaded to the Broad Institute artifactory at https://broadinstitute.jfrog.io/broadinstitute/.
    * Consider updating to the official newer Slick version once they fix issue #2076
    * Related Slick PR: https://github.com/slick/slick/pull/2101
    *
    * Update 2022-03-23: This #2201 PR cherry picks Greg's #2101 PR above and claims to fix the issue:
    * https://github.com/slick/slick/pull/2201
  */
  private val slickV = "3.4.0-M1"
  private val snakeyamlV = "1.33"
  private val sprayJsonV = "1.3.6"
  private val sttpV = "1.7.2"
  private val swaggerParserV = "1.0.56"
  private val swaggerUiV = "4.5.2"
  private val testContainersScalaV = "0.40.10"
  private val tikaV = "2.3.0"
  private val typesafeConfigV = "1.4.2"
  private val workbenchGoogleV = "0.21-5c9c4f6" // via: https://github.com/broadinstitute/workbench-libs/blob/develop/google/CHANGELOG.md
  private val workbenchModelV = "0.15-f9f0d4c" // via: https://github.com/broadinstitute/workbench-libs/blob/develop/model/CHANGELOG.md
  private val workbenchUtilV = "0.6-65bba14" // via: https://github.com/broadinstitute/workbench-libs/blob/develop/util/CHANGELOG.md

  private val slf4jFacadeDependencies = List(
    "org.slf4j" % "slf4j-api" % slf4jV,
    "com.typesafe.scala-logging" %% "scala-logging" % scalaLoggingV,
  )

  private val circeYamlDependency = "io.circe" %% "circe-yaml" % circeYamlV

  private val circeDependencies = List(
    "core",
    "parser",
    "generic",
    "shapes",
    "refined",
    "literal",
    "optics"
  ).map(m => "io.circe" %% s"circe-$m" % circeV) :+ circeYamlDependency :+
  "io.circe" %% "circe-generic-extras" % circeGenericExtrasV :+
  "io.circe" %% "circe-config" % circeConfigV

  private val catsDependencies = List(
    "org.typelevel" %% "cats-core" % catsV,
    "org.typelevel" %% "alleycats-core" % catsV,
    "org.typelevel" %% "mouse" % mouseV,
    "org.typelevel" %% "kittens" % kittensV
  )

  private val http4sDependencies = List(
    "org.http4s" %% "http4s-dsl" % http4sV,
    "org.http4s" %% "http4s-blaze-client" % http4sV,
    "org.http4s" %% "http4s-circe" % http4sV,
  )

  private val googleApiClientDependencies = List(
    // Used by swagger, but only in tests.  This overrides an older 2.1.3 version of jackson-core brought in by
    // these Google dependencies, but which isn't properly evicted by IntelliJ's sbt integration.
    "com.fasterxml.jackson.core" % "jackson-core" % jacksonV,
    // The exclusions prevent guava from colliding at assembly time.
    "com.google.guava" % "guava" % guavaV,
    "com.google.api-client" % "google-api-client-java6" % googleApiClientV
      exclude("com.google.guava", "guava-jdk5"),
    "com.google.api-client" % "google-api-client-jackson2" % googleApiClientV
      exclude("com.google.guava", "guava-jdk5"),
    "com.google.cloud" % "google-cloud-resourcemanager" % googleCloudResourceManagerV,
    /*
    The google-cloud-java dependencies have similar issues with using an older javax.* vs. jakarta.* as guice.
    google-cloud-java is still using javax.annotation and guice is sticking with javax.inject:
     - https://github.com/google/guice/issues/1383
     - https://github.com/googleapis/google-cloud-java/blob/v0.201.0/google-cloud-jar-parent/pom.xml#L131-L136

     Globally use of jakarta instead of javax until Google does themselves.
     The javax.annotation exclusion is below in cromwellExcludeDependencies.
     */
    "jakarta.annotation" % "jakarta.annotation-api" % jakartaAnnotationV,
  )

  val spiDependencies: List[ModuleID] = List(
    "com.iheart" %% "ficus" % ficusV,
  ) ++ googleApiClientDependencies ++ slf4jFacadeDependencies

  val spiUtilDependencies = List(
    "com.iheart" %% "ficus" % ficusV,
    "org.typelevel" %% "cats-effect" % catsEffectV,
  )

  /*
   Due to complications with the public preview Nio filesystem for azure,
  we include this FS implementation locally and include its dependencies
   */
  val azureBlobNioDependencies = List(
    "com.azure" % "azure-core" % "1.40.0",
    "com.azure" % "azure-storage-blob" % "12.23.0-beta.1",
    "com.azure" % "azure-storage-common" % "12.22.0-beta.1",
    "com.azure" % "azure-core-test" % "1.18.0",
    "org.junit.jupiter" % "junit-jupiter-params" % "5.9.3",
    "org.junit.jupiter" % "junit-jupiter-engine" % "5.9.3",
    "org.junit.jupiter" % "junit-jupiter-api" % "5.9.3",
    "io.projectreactor" % "reactor-test" % "3.4.29",
    "cglib" % "cglib-nodep" % "3.2.7",
    "com.azure" % "azure-core-http-okhttp" % "1.11.10",
    "org.mockito" % "mockito-core" % "4.11.0",
    "com.github.sbt" % "junit-interface" % "0.13.2"
  )

  val azureDependencies: List[ModuleID] = List(
    "com.azure" % "azure-identity" % azureIdentitySdkV
      exclude("jakarta.xml.bind", "jakarta.xml.bind-api")
      exclude("jakarta.activation", "jakarta.activation-api")
      exclude("net.minidev", "json-smart"),
    "com.azure" % "azure-identity-extensions" % azureIdentityExtensionsV,
    "com.azure" % "azure-core-management" % azureCoreManagementV,
    "com.fasterxml.jackson.dataformat" % "jackson-dataformat-xml" % jacksonV,
    "com.azure.resourcemanager" % "azure-resourcemanager" % "2.18.0",
    "net.minidev" % "json-smart" % jsonSmartV,
    "com.microsoft.azure" % "applicationinsights-logging-logback" % azureAppInsightsLogbackV,
  ) ++ azureBlobNioDependencies

  val wsmDependencies: List[ModuleID] = List(
    "bio.terra" % "workspace-manager-client" % "0.254.452-SNAPSHOT"
      exclude("com.sun.activation", "jakarta.activation"),
    /*
    1. WSM is looking for the rs-api under javax.*.

    Jersey 3.x switched to jakarta.ws.rs-api 3.x. If one uses jakarta's rs-api, 3.x will automatically evict 2.x.

    However, jakarta's rs-api 2.x provides packages javax.* while 3.x provides jakarta.* instead.
     - https://javadoc.io/doc/jakarta.ws.rs/jakarta.ws.rs-api/2.1.6/javax/ws/rs/package-summary.html
     - https://javadoc.io/doc/jakarta.ws.rs/jakarta.ws.rs-api/3.1.0/jakarta.ws.rs/module-summary.html

    TODO: Perhaps coordinate with the WSM team to use the jakarta 3.x rs-api and jakarta.* instead of javax.*.

    2. Use the exact version of jersey that WSM is using.

    Jersey libraries cannot be mixed and matched as the various modules cannot be mixed and matched.
    For example jersey-client 2.32 is not compatible with jersey-common 2.37.

    If needed one may also explicitly enumerate the list of jersey artifacts and explicitly set the versions similar to
    catsDepeendencies, akkaHttpDependencies, etc.
     - https://broadinstitute.jfrog.io/ui/repos/tree/PomView/libs-snapshot-local/bio/terra/workspace-manager-client/0.254.452-SNAPSHOT/workspace-manager-client-0.254.452-20221114.190249-1.pom
     - https://github.com/eclipse-ee4j/jersey/blob/2.32/core-client/src/main/java/org/glassfish/jersey/client/ClientExecutorProvidersConfigurator.java#L139
     - https://github.com/eclipse-ee4j/jersey/blob/2.37/core-client/src/main/java/org/glassfish/jersey/client/ClientExecutorProvidersConfigurator.java#L136-L137
     */
    "org.glassfish.jersey.inject" % "jersey-hk2" % jerseyV
      exclude("com.sun.activation", "jakarta.activation"),
    "jakarta.activation" % "jakarta.activation-api" % jakartaActivationV,
  )

  val implFtpDependencies = List(
    "commons-net" % "commons-net" % commonNetV,
    "io.github.andrebeat" %% "scala-pool" % scalaPoolV,
    "com.google.guava" % "guava" % guavaV,
    "org.mockftpserver" % "MockFtpServer" % mockFtpServerV % Test
  )

  val implDrsDependencies: List[ModuleID] = List(
    "org.apache.commons" % "commons-lang3" % commonsLang3V,
    "com.google.cloud" % "google-cloud-storage" % googleCloudStorageV,
    "com.google.oauth-client" % "google-oauth-client" % googleOauthClientV
  ) ++ circeDependencies ++ catsDependencies ++ azureDependencies

  // Internal collections of dependencies

  private val betterFilesDependencies = List(
    "com.github.pathikrit" %% "better-files" % betterFilesV
  )

  private val configDependencies = List(
    "com.typesafe" % "config" % typesafeConfigV,
    "com.iheart" %% "ficus" % ficusV
  )

  /*
  Adds a variety of logging libraries required for actual logging. However, some of these aren't always required.

  Ex: If one isn't using akka & slf4j, then 'akka-slf4j' isn't required. However, for now, all executables are using
  akka & slf4j... so leaving it.

  Similarly, not _all_ executables/logback.xml configs will need logback-access, raven-logback, janino, etc.
  Still, leaving them as dependencies for simplicity's sake.
   */
  private val slf4jBindingDependencies = List(
    // http://logback.qos.ch/dependencies.html
    "ch.qos.logback" % "logback-access" % logbackV,
    "ch.qos.logback" % "logback-classic" % logbackV,
    "ch.qos.logback" % "logback-core" % logbackV,
    "com.typesafe.akka" %% "akka-slf4j" % akkaV,
    "io.sentry" % "sentry-logback" % sentryLogbackV,
    "org.codehaus.janino" % "janino" % janinoV,
    // Replace all log4j usage with slf4j
    // https://www.slf4j.org/legacy.html#log4j-over-slf4j
    "org.slf4j" % "log4j-over-slf4j" % slf4jV,
    // Replace all commons-logging usage with slf4j
    // https://www.slf4j.org/legacy.html#jcl-over-slf4j
    "org.slf4j" % "jcl-over-slf4j" % slf4jV,
    // Enable runtime replacing of java.util.logging usage with slf4j
    // https://www.slf4j.org/legacy.html#jul-to-slf4j
    "org.slf4j" % "jul-to-slf4j" % slf4jV,
  ) ++ slf4jFacadeDependencies

  private val slickDependencies = List(
    "com.typesafe.slick" %% "slick" % slickV,
    "com.typesafe.slick" %% "slick-hikaricp" % slickV,
    "com.rms.miu" %% "slick-cats" % slickCatsV
  )

  private val liquibaseDependencies = List(
    // The XML bind API replacement below may be removed when this ticket is addressed:
    // https://github.com/liquibase/liquibase/issues/2991
    "org.liquibase" % "liquibase-core" % liquibaseV
      exclude("javax.xml.bind", "jaxb-api"),
    "jakarta.xml.bind" % "jakarta.xml.bind-api" % jakartaXmlBindApiV,
  )

  private val akkaDependencies = List(
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-testkit" % akkaV % Test,
  )

  private val akkaStreamDependencies = List(
    "com.typesafe.akka" %% "akka-stream" % akkaV,
    "com.typesafe.akka" %% "akka-stream-testkit" % akkaV % Test,
  ) ++ akkaDependencies

  private val akkaHttpDependencies = List(
    "com.typesafe.akka" %% "akka-http" % akkaHttpV,
    "com.typesafe.akka" %% "akka-http-testkit" % akkaHttpV % Test,
    // WOM internally embeds spray-json. Leave this import here until WOM externalizes the json library choice like
    // other libraries do. See akka-http, elastic4s, etc.
    "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpV,
  ) ++ akkaStreamDependencies

  private val akkaHttpCirceIntegrationDependency = List(
    "de.heikoseeberger" %% "akka-http-circe" % akkaHttpCirceIntegrationV
  )

  private val swaggerUiDependencies = List(
    "org.webjars" % "swagger-ui" % swaggerUiV,
    "io.swagger" % "swagger-parser" % swaggerParserV % Test,
    "org.yaml" % "snakeyaml" % snakeyamlV % Test
  )

  private val googleGenomicsV2Alpha1Dependency = List(
    "com.google.apis" % "google-api-services-genomics" % googleGenomicsServicesV2Alpha1ApiV
      exclude("com.google.guava", "guava-jdk5")
  )

  private val googleLifeSciencesV2BetaDependency = List(
    "com.google.apis" % "google-api-services-lifesciences" % googleLifeSciencesServicesV2BetaApiV
      exclude("com.google.guava", "guava-jdk5")
  )

  private val googleBatchv1Dependency = List(
    "com.google.cloud" % "google-cloud-batch" % googleCloudBatchV1,
    "com.google.api.grpc" % "proto-google-cloud-batch-v1" % googleCloudBatchV1,
    "com.google.api.grpc" % "proto-google-cloud-resourcemanager-v3" % "1.17.0"
  )

  /*
  Used instead of `"org.lerch" % "s3fs" % s3fsV exclude("org.slf4j", "jcl-over-slf4j")`
  org.lerch:s3fs:1.0.1 depends on a preview release of software.amazon.awssdk:s3.

  Instead the code has been re-forked into this repo, just like many of the other FileSystemProvider extensions.
   */
  private val s3fsDependencies = List(
    "com.google.code.findbugs" % "jsr305" % jsr305V,
    "com.google.guava" % "guava" % guavaV,
    "org.apache.tika" % "tika-core" % tikaV,
    "software.amazon.awssdk" % "s3" % awsSdkV,
  ) ++ slf4jBindingDependencies

  private val awsCloudDependencies = List(
    "com.fasterxml.jackson.core" % "jackson-annotations" % jacksonV,
  ) ++ s3fsDependencies ++ List(
    "batch",
    "core",
    "cloudwatchlogs",
    "ecr",
    "ecrpublic",
    "s3",
    "sts",
  ).map(artifactName => "software.amazon.awssdk" % artifactName % awsSdkV)

  private val googleCloudDependencies = List(
    "io.grpc" % "grpc-core" % grpcV,
    "com.google.guava" % "guava" % guavaV,
    /*
    The google-cloud-nio has the same problems with an ancient inject as guice:
     - https://github.com/google/guice/issues/1383
     - https://github.com/googleapis/java-storage-nio/blob/v0.124.20/google-cloud-nio/pom.xml#L49-L53

     Force use of jakarta instead of javax until Google does themselves.
     */
    "com.google.cloud" % "google-cloud-nio" % googleCloudNioV
      exclude("com.google.api.grpc", "grpc-google-common-protos")
      exclude("com.google.cloud.datastore", "datastore-v1-protos")
      exclude("javax.inject", "javax.inject")
      exclude("org.apache.httpcomponents", "httpclient"),
    "org.broadinstitute.dsde.workbench" %% "workbench-google" % workbenchGoogleV
      exclude("com.google.apis", "google-api-services-genomics"),
    "org.apache.httpcomponents.client5" % "httpclient5" % apacheHttpClient5V,
  "com.google.apis" % "google-api-services-cloudkms" % googleCloudKmsV
      exclude("com.google.guava", "guava-jdk5"),
    "org.glassfish.hk2.external" % "jakarta.inject" % jakartaInjectV,
  ) ++ googleGenomicsV2Alpha1Dependency ++ googleLifeSciencesV2BetaDependency ++ googleBatchv1Dependency

  private val dbmsDependencies = List(
    "org.hsqldb" % "hsqldb" % hsqldbV,
    "org.mariadb.jdbc" % "mariadb-java-client" % mariadbV,
    "mysql" % "mysql-connector-java" % mysqlV,
    "org.postgresql" % "postgresql" % postgresV
  )

  private val refinedTypeDependenciesList = List(
    "eu.timepit" %% "refined" % refinedV
  )

  // Sub-project dependencies, added in addition to any dependencies inherited from .dependsOn().

  val commonDependencies: List[ModuleID] = List(
    "org.typelevel" %% "cats-effect" % catsEffectV,
    "org.apache.commons" % "commons-lang3" % commonsLang3V,
    "org.apache.commons" % "commons-text" % commonsTextV,
    "com.lihaoyi" %% "pprint" % pprintV,
  ) ++ catsDependencies ++ configDependencies ++ slf4jFacadeDependencies ++ refinedTypeDependenciesList

  val cloudSupportDependencies: List[ModuleID] = googleApiClientDependencies ++ googleCloudDependencies ++ betterFilesDependencies ++ awsCloudDependencies ++ azureDependencies

  val databaseSqlDependencies: List[ModuleID] = List(
    "commons-io" % "commons-io" % commonsIoV,
  ) ++ configDependencies ++ catsDependencies ++ slickDependencies ++ dbmsDependencies ++ refinedTypeDependenciesList

  val statsDDependencies = List(
    "nl.grons" %% "metrics4-scala" % metrics4ScalaV,
    "com.readytalk" % "metrics3-statsd" % metrics3StatsdV
  )

  val stackdriverDependencies = List(
    "com.google.cloud" % "google-cloud-monitoring" % googleCloudMonitoringV
  )

  /*
  Generators are eventually coming to ScalaTest. Someday...
    - https://youtu.be/lKtg-CDVDsI?t=562

  For now use scalatestplus' scalacheck wrapper.

  Tests that insist on using PropertyGenerators should actually use ScalaTest's wrapper. ScalaCheck tests no longer
  run by default. See Testing.scala where only `ScalaTest` is specified in the `testFrameworks`.

  See also (may be out of date):
    - https://github.com/scalatest/scalatest/issues/1735
    - https://www.scalatest.org/user_guide/generator_driven_property_checks
    - https://www.scalatest.org/user_guide/writing_scalacheck_style_properties
   */
  private val scalacheckBaseV = "1.17"
  private val scalacheckDependencies = List(
    "org.scalatestplus" %% s"scalacheck-${scalacheckBaseV.replace(".", "-")}" % scalatestScalacheckV % Test,
  )

  /*
  Note: `junitDependencies` only adds the dependency for JUnit tests to compile.

  To actually _run_ the tests via SBT one would need the SBT to JUnit interface:
    - https://github.com/sbt/junit-interface/

  However, as of Aug 2021 there is only one S3 Java file using JUnit, and that code was copy-pasted from an
  external GitHub repo. See `s3fsDependencies` for more information.

  Also as of Aug 2021 Testing.scala only looks for and runs ScalaTest during regular testing.
   */
  private val junitDependencies = List(
    "junit" % "junit" % junitV % Test
  )

  private val testDatabaseDependencies =
    List("scalatest", "mysql", "mariadb", "postgresql")
      .map(name => "com.dimafeng" %% s"testcontainers-scala-$name" % testContainersScalaV % Test)

  val blobFileSystemDependencies: List[ModuleID] = azureDependencies ++ wsmDependencies ++ akkaHttpDependencies

  val s3FileSystemDependencies: List[ModuleID] = junitDependencies

  val gcsFileSystemDependencies: List[ModuleID] = akkaHttpDependencies

  val httpFileSystemDependencies: List[ModuleID] = akkaHttpDependencies

  val womDependencies: List[ModuleID] = List(
    "com.typesafe.scala-logging" %% "scala-logging" % scalaLoggingV,
    "io.spray" %% "spray-json" % sprayJsonV,
    "org.typelevel" %% "simulacrum" % simulacrumV,
    "commons-codec" % "commons-codec" % commonsCodecV
  ) ++ scalacheckDependencies ++ circeDependencies ++ refinedTypeDependenciesList

  val wdlDependencies: List[ModuleID] = List(
    "commons-io" % "commons-io" % commonsIoV,
    "org.scala-graph" %% "graph-core" % scalaGraphV,
    "com.chuusai" %% "shapeless" % shapelessV,
    "com.google.re2j" % "re2j" % re2jV,
  ) ++ betterFilesDependencies

  val languageFactoryDependencies = List(
    "com.softwaremill.sttp" %% "core" % sttpV,
    "com.softwaremill.sttp" %% "async-http-client-backend-cats" % sttpV
      exclude("com.sun.activation", "javax.activation"),
    "jakarta.activation" % "jakarta.activation-api" % jakartaActivationV,
  )

  val mockServerDependencies = List(
    "org.mock-server" % "mockserver-netty" % mockserverNettyV % Test
  )

  val womtoolDependencies: List[ModuleID] = catsDependencies ++ slf4jBindingDependencies

  val coreDependencies: List[ModuleID] = List(
    "com.google.auth" % "google-auth-library-oauth2-http" % googleOauth2V,
    "com.chuusai" %% "shapeless" % shapelessV,
    // NOTE: See scalameter comment under engineDependencies
    "com.storm-enroute" %% "scalameter" % scalameterV % Test
      exclude("com.fasterxml.jackson.module", "jackson-module-scala_2.13")
      exclude("org.scala-lang.modules", "scala-xml_2.13"),
    "com.github.scopt" %% "scopt" % scoptV,
  ) ++ akkaStreamDependencies ++ configDependencies ++ catsDependencies ++ circeDependencies ++
    googleApiClientDependencies ++ statsDDependencies ++ betterFilesDependencies ++
    // TODO: We're not using the "F" in slf4j. Core only supports logback, specifically the WorkflowLogger.
    slf4jBindingDependencies ++ stackdriverDependencies

  val databaseMigrationDependencies: List[ModuleID] = liquibaseDependencies ++ dbmsDependencies

  val dockerHashingDependencies: List[ModuleID] = http4sDependencies ++ circeDependencies

  val cromwellApiClientDependencies: List[ModuleID] = List(
    "org.typelevel" %% "cats-effect" % catsEffectV,
    "co.fs2" %% "fs2-io" % fs2V % Test,
  ) ++ akkaHttpDependencies ++ betterFilesDependencies ++ catsDependencies

  val centaurDependencies: List[ModuleID] = List(
    "org.apache.commons" % "commons-math3" % commonsMathV,
    "com.github.kxbmap" %% "configs" % configsV,
    "com.google.cloud" % "google-cloud-bigquery" % googleCloudBigQueryV % IntegrationTest,
    "org.gnieh" %% "diffson-spray-json" % diffsonSprayJsonV
  ) ++ circeDependencies ++ slf4jBindingDependencies ++ cloudSupportDependencies ++ http4sDependencies

  val engineDependencies: List[ModuleID] = List(
    "commons-codec" % "commons-codec" % commonsCodecV,
    "commons-io" % "commons-io" % commonsIoV,
    /*
    Maybe ScalaMeter should be used, but is anyone?

    For now keep its dependencies from breaking jackson for other libraries. If someone wants to use it they can
    re-fight with dependency-hell at that point.

    Avoid:
    "com.fasterxml.jackson.databind.JsonMappingException: Scala module 2.11.3 requires Jackson Databind
    version >= 2.11.0 and < 2.12.0":
     - https://scalameter.github.io/home/gettingstarted/0.7/sbt/index.html
     - https://github.com/FasterXML/jackson-module-scala/blob/jackson-module-scala-2.11.3/src/main/scala/com/fasterxml/jackson/module/scala/JacksonModule.scala#L53-L62
     */
    "com.storm-enroute" %% "scalameter" % scalameterV
      exclude("com.fasterxml.jackson.core", "jackson-databind")
      exclude("com.fasterxml.jackson.module", "jackson-module-scala")
      exclude("com.fasterxml.jackson.module", "jackson-module-scala_2.13")
      exclude("org.scala-tools.testing", "test-interface")
      exclude("org.scala-lang.modules", "scala-xml_2.13"),
    "com.fasterxml.jackson.core" % "jackson-databind" % jacksonV,
    "io.github.andrebeat" %% "scala-pool" % scalaPoolV
  ) ++ swaggerUiDependencies ++ akkaHttpDependencies ++ akkaHttpCirceIntegrationDependency ++ circeDependencies ++
    testDatabaseDependencies

  val servicesDependencies: List[ModuleID] = List(
    "com.google.api" % "gax-grpc" % googleGaxGrpcV,
    "org.apache.commons" % "commons-csv" % commonsCsvV,
    "bio.terra" % "bard-client-resttemplate-javax" % bardClientV
      exclude("org.springframework", "spring-aop")
      exclude("org.springframework", "spring-jcl"),
    "org.apache.httpcomponents.client5" % "httpclient5" % apacheHttpClient5V // Needed for rest-template connection pooling

  ) ++ testDatabaseDependencies ++ akkaHttpDependencies ++ mockServerDependencies

  val serverDependencies: List[ModuleID] = slf4jBindingDependencies

  val cromiamDependencies: List[ModuleID] = List(
    /*
    sttp 1.x was last released in 2019
    See above comment regarding "cats-effect, fs2, http4s, and sttp" all needing to update together.
    For now, replace sttp 1.x's com.sun.activation usage with the jakarta version.

    NOTE when upgrading: sttp 3.x no longer requires an async-http-client-backend-future so jakarta.activation can
    probably be removed from the dependencies:
     - https://sttp.softwaremill.com/en/v3/backends/future.html#using-async-http-client
     - https://sttp.softwaremill.com/en/v2/backends/future.html#using-async-http-client
     - https://sttp.softwaremill.com/en/v1/backends/asynchttpclient.html
     */
    "com.softwaremill.sttp" %% "core" % sttpV,
    "com.softwaremill.sttp" %% "async-http-client-backend-future" % sttpV
      exclude("com.sun.activation", "javax.activation"),
    "com.typesafe.scala-logging" %% "scala-logging" % scalaLoggingV,
    "org.broadinstitute.dsde.workbench" %% "workbench-model" % workbenchModelV,
    "org.broadinstitute.dsde.workbench" %% "workbench-util" % workbenchUtilV,
    "jakarta.activation" % "jakarta.activation-api" % jakartaActivationV,
  ) ++ akkaHttpDependencies ++ swaggerUiDependencies ++ slf4jBindingDependencies

  val wes2cromwellDependencies: List[ModuleID] = coreDependencies ++ akkaHttpDependencies

  val backendDependencies: List[ModuleID] = List(
    "co.fs2" %% "fs2-io" % fs2V
  ) ++ scalacheckDependencies

  val tesBackendDependencies: List[ModuleID] = akkaHttpDependencies

  val sfsBackendDependencies = List (
    "org.lz4" % "lz4-java" % lz4JavaV
  )
  val scalaTest = "org.scalatest" %% "scalatest" % scalatestV

  val testDependencies: List[ModuleID] = List(
    "org.scalatest" %% "scalatest" % scalatestV,
    // Use mockito Java DSL directly instead of the numerous and often hard to keep updated Scala DSLs.
    // See also scaladoc in common.mock.MockSugar and that trait's various usages.
    "org.mockito" % "mockito-core" % mockitoV
  ) ++ slf4jBindingDependencies // During testing, add an slf4j binding for _all_ libraries.

  val kindProjectorPlugin = "org.typelevel" % "kind-projector" % kindProjectorV cross CrossVersion.full

  // Version of the swagger UI to write into config files
  val swaggerUiVersion: String = swaggerUiV

  val drsLocalizerDependencies: List[ModuleID] = List(
    "com.google.auth" % "google-auth-library-oauth2-http" % googleOauth2V,
    "com.google.cloud" % "google-cloud-storage" % googleCloudStorageV,
    "org.typelevel" %% "cats-effect" % catsEffectV,
    "com.iheart" %% "ficus" % ficusV,
    "com.softwaremill.sttp" %% "circe" % sttpV,
    "com.github.scopt" %% "scopt" % scoptV,
    "org.apache.commons" % "commons-csv" % commonsCsvV,
    "io.spray" %% "spray-json" % sprayJsonV,
  ) ++ circeDependencies ++ catsDependencies ++ slf4jBindingDependencies ++ languageFactoryDependencies ++ azureDependencies

  val allProjectDependencies: List[ModuleID] =
    backendDependencies ++
      centaurDependencies ++
      cloudSupportDependencies ++
      commonDependencies ++
      coreDependencies ++
      cromiamDependencies ++
      cromwellApiClientDependencies ++
      databaseMigrationDependencies ++
      databaseSqlDependencies ++
      mockServerDependencies ++
      drsLocalizerDependencies ++
      engineDependencies ++
      gcsFileSystemDependencies ++
      httpFileSystemDependencies ++
      implDrsDependencies ++
      implFtpDependencies ++
      languageFactoryDependencies ++
      serverDependencies ++
      sfsBackendDependencies ++
      spiDependencies ++
      spiUtilDependencies ++
      tesBackendDependencies ++
      wdlDependencies ++
      wes2cromwellDependencies ++
      womDependencies ++
      womtoolDependencies

  /*
  If you see warnings from SBT about evictions, insert a specific dependency version into this list.

  Do not know a good way to check when these are out of date as `sbt dependencyUpdates` does not
  report on dependency overrides.

  Any dependencies that are removed may be also removed from this list.
  However, be careful about downgrading any of these dependencies.
  Older versions have known vulnerabilities, ex: CVE-2017-7525

  === SECURITY UPGRADES ===

  When upgrading dependencies to fix security issues, it is preferable to start with upgrading the
  library that brings it in. Only fall back to overriding here when the latest library version still
  has a vulnerable version of the dependency, or a major version upgrade is required and infeasible.
  This algorithm makes it simpler to upgrade libraries in the future, because we don't have to
  remember to remove the override.
   */

  val googleHttpClientDependencies = List(
    /*
    Move the google-http-client versions past https://github.com/googleapis/google-http-java-client/issues/606
    This created a situation where com/google/api/client/http/apache/ApacheHttpTransport.class was in *both*
    transitive dependencies causing an assembly merge conflict.

    At the time of this comment older versions are being pulled in via
    https://mvnrepository.com/artifact/com.google.api-client/google-api-client/1.28.0
     */
    "com.google.http-client" % "google-http-client-apache" % googleHttpClientApacheV,
    "com.google.http-client" % "google-http-client" % googleHttpClientV,
  )

  val nettyDependencyOverrides: List[ModuleID] = List(
    "buffer",
    "codec",
    "codec-dns",
    "codec-http",
    "codec-http2",
    "codec-socks",
    "common",
    "handler-proxy",
    "resolver",
    "resolver-dns",
    "transport",
    "transport-native-epoll",
    "transport-native-unix-common",
  ).map(m => "io.netty" % s"netty-$m" % nettyV)

  val rdf4jDependencyOverrides: List[ModuleID] = List(
    /*
    Yes. All of these are required to lock in the rdf4j version.

    Feel free to update versions but do not remove these overrides unless and until an updated
    owl-api is no longer pulling in vulnerable rdf4j dependencies.

    https://cve.mitre.org/cgi-bin/cvename.cgi?name=CVE-2018-1000644

    See comment mentioning "OSGI" further above for more info on the bundling of dependencies.
     */
    "model",
    "rio-api",
    "rio-binary",
    "rio-datatypes",
    "rio-jsonld",
    "rio-languages",
    "rio-n3",
    "rio-nquads",
    "rio-ntriples",
    "rio-rdfjson",
    "rio-rdfxml",
    "rio-trig",
    "rio-trix",
    "rio-turtle",
    "util",
  ).map(m => "org.eclipse.rdf4j" % s"rdf4j-$m" % rdf4jV)

  // Some libraries are importing older version of these dependencies, causing conflicts. Hence the need to override them.
  val grpcDependencyOverrides: List[ModuleID] = List(
    "alts",
    "auth",
    "context",
    "core",
    "grpclb",
    "netty-shaded",
    "protobuf-lite",
    "protobuf",
    "stub",
  ).map(m => "io.grpc" % s"grpc-$m" % grpcV)

  /*
  Ensure we're using the latest to avoid a shading bug in earlier versions of scala-collection-compat.
  https://github.com/scala/scala-collection-compat/issues/426
   */
  private val scalaCollectionCompatOverrides = List(
    "org.scala-lang.modules" %% "scala-collection-compat" % scalaCollectionCompatV,
  )

  private val asyncHttpClientOverrides = List(
    "org.asynchttpclient" % "async-http-client" % "2.10.5",
  )


  private val nimbusdsOverrides = List(
    "com.nimbusds" % "nimbus-jose-jwt" % "9.23",
  )

  private val bouncyCastleOverrides = List(
    "org.bouncycastle" % "bcprov-jdk15on" % "1.70",
  )

  private val protobufJavaOverrides = List(
    "com.google.protobuf" % "protobuf-java" % "3.21.2",
  )

  /*
  If we use a version in one of our projects, that's the one we want all the libraries to use
  ...plus other groups of transitive dependencies shared across multiple projects
   */
  val cromwellDependencyOverrides: List[ModuleID] =
    allProjectDependencies ++
      googleHttpClientDependencies ++
      nettyDependencyOverrides ++
      rdf4jDependencyOverrides ++
      grpcDependencyOverrides ++
      scalaCollectionCompatOverrides ++
      asyncHttpClientOverrides ++
      nimbusdsOverrides ++
      bouncyCastleOverrides ++
      protobufJavaOverrides

  /*
  Libraries that should be globally excluded.
   */
  val cromwellExcludeDependencies: List[ExclusionRule] = List(
    // Replaced with jcl-over-slf4j
    ExclusionRule("commons-logging", "commons-logging"),
    /*
    The google-cloud-java dependencies have similar issues with using an older javax.* vs. jakarta.* as guice.
    google-cloud-java is still using javax.annotation and guice is sticking with javax.inject:
     - https://github.com/google/guice/issues/1383
     - https://github.com/googleapis/google-cloud-java/blob/v0.201.0/google-cloud-jar-parent/pom.xml#L131-L136

     Globally use of jakarta instead of javax until Google does themselves.
     The jakarta.annotation inclusion is above in googleApiClientDependencies.
     */
    ExclusionRule("javax.annotation", "javax.annotation-api"),
    ExclusionRule("javax.activation"),
  )

  val http4sDsl = "org.http4s" %% "http4s-dsl" % http4sV
  val http4sEmberClient = "org.http4s" %% "http4s-ember-client" % http4sV
  val http4sEmberServer = "org.http4s" %% "http4s-ember-server" % http4sV
  val http4sCirce = "org.http4s" %% "http4s-circe" % http4sV
  val pact4sScalaTest = "io.github.jbwheatley" %% "pact4s-scalatest" % pact4sV % Test
  val pact4sCirce = "io.github.jbwheatley" %% "pact4s-circe" % pact4sV
  val pact4sSpray = "io.github.jbwheatley" %% "pact4s-spray-json" % pact4sV

  val pact4sDependencies = Seq(
    pact4sScalaTest,
    pact4sCirce,
    pact4sSpray,
    http4sEmberClient,
    http4sDsl,
    http4sEmberServer,
    http4sCirce,
    scalaTest,
  ) ++ akkaDependencies
}
