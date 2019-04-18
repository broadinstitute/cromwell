import sbt._

object Dependencies {
  private val akkaHttpCirceIntegrationV = "1.24.3"
  private val akkaHttpV = "10.1.7"
  private val akkaV = "2.5.19"
  private val aliyunBcsV = "6.0.6"
  private val aliyunCoreV = "4.3.2"
  private val aliyunOssV = "3.4.0"
  private val ammoniteOpsV = "1.6.3"
  private val apacheCommonNetV = "3.6"
  private val apacheHttpClientV = "4.5.7"
  private val awsSdkV = "2.3.9"
  private val betterFilesV = "2.17.1"
  private val catsEffectV = "1.2.0"
  private val catsV = "1.5.0"
  private val circeOpticsV = "0.11.0"
  private val circeV = "0.11.1"
  private val circeYamlV = "0.9.0"
  private val commonsCodecV = "1.11"
  private val commonsIoV = "2.6"
  private val commonsLang3V = "3.8.1"
  private val commonsTextV = "1.6"
  private val configsV = "0.4.4"
  private val delightRhinoSandboxV = "0.0.10"
  private val ficusV = "1.4.4"
  private val fs2V = "1.0.3"
  private val googleApiClientV = "1.28.0"
  private val googleCloudCoreV = "1.61.0"
  private val googleCloudKmsV = "v1-rev63-1.25.0"
  private val googleCloudNioV = "0.61.0-alpha"
  private val googleGenomicsServicesV1ApiV = "v1alpha2-rev495-1.23.0"
  private val googleGenomicsServicesV2ApiV = "v2alpha1-rev31-1.25.0"
  private val googleOauth2V = "0.13.0"
  private val grpcV = "1.18.0"
  private val guavaV = "27.0.1-jre"
  private val heterodonV = "1.0.0-beta3"
  private val hsqldbV = "2.4.1"
  private val http4sVersion = "0.20.0-M5"
  private val jacksonV = "2.9.8"
  private val janinoV = "3.0.12"
  private val javaxActivationV = "1.2.0"
  private val jaxbV = "2.3.2"
  private val kindProjectorV = "0.9.9"
  private val kittensV = "1.2.0"
  private val liquibaseSlf4jV = "2.0.0"
  private val liquibaseV = "3.5.5" // https://github.com/broadinstitute/cromwell/issues/4618
  private val logbackV = "1.2.3"
  private val metrics3ScalaV = "3.5.10" // https://github.com/erikvanoosten/metrics-scala/tree/f733e26#download-4x
  private val metrics3StatsdV = "4.2.0"
  private val mockFtpServerV = "2.7.1"
  private val mockserverNettyV = "5.5.1"
  private val mouseV = "0.20"
  private val mysqlV = "8.0.15"
  private val nettyV = "4.1.33.Final"
  private val owlApiV = "5.1.9"
  private val paradiseV = "2.1.1"
  private val pegdownV = "1.6.0"
  private val rdf4jV = "2.4.2"
  private val refinedV = "0.9.4"
  private val rhinoV = "1.7.10"
  private val scalaGraphV = "1.12.5"
  private val scalaLoggingV = "3.9.2"
  private val scalaPoolV = "0.4.1"
  private val scalacheckV = "1.14.0"
  private val scalacticV = "3.0.5"
  private val scalameterV = "0.10.1"
  private val scalamockV = "4.1.0"
  private val scalatestV = "3.0.5"
  private val scalazV = "7.2.27"
  private val scoptV = "3.7.1"
  private val sentryLogbackV = "1.7.17"
  private val shapelessV = "2.3.3"
  private val simulacrumV = "0.15.0"
  private val slf4jV = "1.7.25"
  private val slickCatsV = "0.9.0"
  private val slickV = "3.2.3"
  private val snakeyamlV = "1.23"
  private val specs2MockV = "4.4.1"
  private val sprayJsonV = "1.3.5"
  private val sttpV = "1.5.8"
  private val swaggerParserV = "1.0.41"
  private val swaggerUiV = "3.2.2"
  private val tikaV = "1.20"
  private val typesafeConfigV = "1.3.3"
  private val workbenchGoogleV = "0.15-2fc79a3"
  private val workbenchModelV = "0.10-6800f3a"
  private val workbenchUtilV = "0.3-f3ce961"

  private val circeYamlDependency = "io.circe" %% "circe-yaml" % circeYamlV

  private val circeDependencies = List(
    "core",
    "parser",
    "generic",
    "generic-extras",
    "java8",
    "shapes",
    "refined",
    "literal"
  ).map(m => "io.circe" %% s"circe-$m" % circeV) :+ circeYamlDependency

  private val catsDependencies = List(
    "org.typelevel" %% "cats-core" % catsV,
    "org.typelevel" %% "alleycats-core" % catsV,
    "org.typelevel" %% "mouse" % mouseV,
    "org.typelevel" %% "kittens" % kittensV
  )

  private val http4sDependencies = List(
    "org.http4s" %% "http4s-dsl" % http4sVersion,
    "org.http4s" %% "http4s-blaze-client" % http4sVersion,
    "org.http4s" %% "http4s-circe" % http4sVersion
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
      exclude("com.google.guava", "guava-jdk5")
  )

  val spiDependencies = List(
    "com.iheart" %% "ficus" % ficusV,
    "org.slf4j" % "slf4j-api" % slf4jV,
    "com.typesafe.scala-logging" %% "scala-logging" % scalaLoggingV
  ) ++ googleApiClientDependencies

  val spiUtilDependencies = List(
    "com.iheart" %% "ficus" % ficusV,
    "org.typelevel" %% "cats-effect" % catsEffectV,
  )

  val implFtpDependencies = List(
    "commons-net" % "commons-net" % apacheCommonNetV,
    "io.github.andrebeat" %% "scala-pool" % scalaPoolV,
    "com.google.guava" % "guava" % guavaV,
    "org.scalamock" %% "scalamock" % scalamockV % Test,
    "org.mockftpserver" % "MockFtpServer" % mockFtpServerV % Test
  )

  val implDrsDependencies = List(
    "org.apache.commons" % "commons-lang3" % commonsLang3V,
    "com.google.cloud" % "google-cloud-storage" % googleCloudCoreV,
    "com.google.oauth-client" % "google-oauth-client" % googleApiClientV
  ) ++ circeDependencies ++ catsDependencies

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
    "org.slf4j" % "log4j-over-slf4j" % slf4jV
  )

  private val slickDependencies = List(
    "com.typesafe.slick" %% "slick" % slickV,
    "com.typesafe.slick" %% "slick-hikaricp" % slickV,
    "com.rms.miu" %% "slick-cats" % slickCatsV
  )

  private val liquibaseDependencies = List(
    "org.liquibase" % "liquibase-core" % liquibaseV,
    // This is to stop liquibase from being so noisy by default
    // See: http://stackoverflow.com/questions/20880783/how-to-get-liquibase-to-log-using-slf4j
    "com.mattbertolini" % "liquibase-slf4j" % liquibaseSlf4jV
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

  // The v1 dependency has been cloned in the broad artifactory so that we can have the 2 versions co-exist in the same jar
  private val googleGenomicsV1Dependency = List(
    "org.broadinstitute" % "cromwell-google-api-services-genomics" % googleGenomicsServicesV1ApiV
      exclude("com.google.guava", "guava-jdk5")
  )

  private val googleGenomicsV2Dependency = List(
    "com.google.apis" % "google-api-services-genomics" % googleGenomicsServicesV2ApiV
      exclude("com.google.guava", "guava-jdk5")
  )

  /*
  Used instead of `"org.lerch" % "s3fs" % s3fsV exclude("org.slf4j", "jcl-over-slf4j")`
  org.lerch:s3fs:1.0.1 depends on a preview release of software.amazon.awssdk:s3.

  Instead the code has been re-forked into this repo, just like many of the other FileSystemProvider extensions.
   */
  private val s3fsDependencies = List(
    "com.google.code.findbugs" % "jsr305" % "3.0.2",
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
    "s3",
    "sts"
  ).map(artifactName => "software.amazon.awssdk" % artifactName % awsSdkV)

  private val googleCloudDependencies = List(
    "io.grpc" % "grpc-core" % grpcV,
    "com.google.guava" % "guava" % guavaV,
    "com.google.cloud" % "google-cloud-nio" % googleCloudNioV
      exclude("com.google.api.grpc", "grpc-google-common-protos")
      exclude("com.google.cloud.datastore", "datastore-v1-protos")
      exclude("org.apache.httpcomponents", "httpclient"),
    "org.broadinstitute.dsde.workbench" %% "workbench-google" % workbenchGoogleV
      exclude("com.google.apis", "google-api-services-genomics"),
    "org.apache.httpcomponents" % "httpclient" % apacheHttpClientV,
    "com.google.apis" % "google-api-services-cloudkms" % googleCloudKmsV
      exclude("com.google.guava", "guava-jdk5")
  ) ++ googleGenomicsV1Dependency ++ googleGenomicsV2Dependency

  private val aliyunOssDependencies = List(
    "com.aliyun.oss" % "aliyun-sdk-oss" % aliyunOssV
      // stax is included twice by oss 3.1.0 and cause assembly merge conflicts via stax vs. javax.xml.stream
      exclude("stax", "stax-api")
      // Exclude jersey-json until aliyun-sdk-oss >3.4.0 is published
      // https://github.com/aliyun/aliyun-oss-java-sdk/pull/149
      exclude("com.sun.jersey", "jersey-json")
      // jaxb-api and jaxb-core and included in jaxb-impl as of 2.3.1
      // https://github.com/eclipse-ee4j/jaxb-ri/issues/1168
      exclude("javax.xml.bind", "jaxb-api")
      exclude("com.sun.xml.bind", "jaxb-core")
      // javax.activation:activation has been replaced. https://stackoverflow.com/a/46493809
      // The old version was causing an assembly merge conflict.
      exclude("javax.activation", "activation"),
    "com.sun.activation" % "javax.activation" % javaxActivationV,
    "com.sun.xml.bind" % "jaxb-impl" % jaxbV,
    "org.glassfish.jaxb" % "jaxb-runtime" % jaxbV
      // already included in com.sun.activation
      exclude("jakarta.activation", "jakarta.activation-api"),
  )

  private val aliyunBatchComputeDependencies = List(
    "com.aliyun" % "aliyun-java-sdk-batchcompute" % aliyunBcsV,
    "com.aliyun" % "aliyun-java-sdk-core" % aliyunCoreV
      // jaxb-api and jaxb-core and included in jaxb-impl as of 2.3.1
      // https://github.com/eclipse-ee4j/jaxb-ri/issues/1168
      exclude("javax.xml.bind", "jaxb-api")
      exclude("com.sun.xml.bind", "jaxb-core")
      // javax.activation:activation has been replaced. https://stackoverflow.com/a/46493809
      // The old version was causing an assembly merge conflict.
      exclude("javax.activation", "activation"),
    "com.sun.activation" % "javax.activation" % javaxActivationV,
    "com.sun.xml.bind" % "jaxb-impl" % jaxbV,
    "org.glassfish.jaxb" % "jaxb-runtime" % jaxbV
      // already included in com.sun.activation
      exclude("jakarta.activation", "jakarta.activation-api"),
  )

  private val dbmsDependencies = List(
    "org.hsqldb" % "hsqldb" % hsqldbV,
    "mysql" % "mysql-connector-java" % mysqlV
  )

  private val refinedTypeDependenciesList = List(
    "eu.timepit" %% "refined" % refinedV
  )

  // Sub-project dependencies, added in addition to any dependencies inherited from .dependsOn().

  val commonDependencies = List(
    "org.slf4j" % "slf4j-api" % slf4jV,
    "org.typelevel" %% "cats-effect" % catsEffectV,
    "org.apache.commons" % "commons-lang3" % commonsLang3V,
    "org.apache.commons" % "commons-text" % commonsTextV,
    "com.typesafe.scala-logging" %% "scala-logging" % scalaLoggingV,
  ) ++ catsDependencies ++ configDependencies

  val cloudSupportDependencies = googleApiClientDependencies ++ googleCloudDependencies ++ betterFilesDependencies ++ awsCloudDependencies

  val databaseSqlDependencies = configDependencies ++ catsDependencies ++ slickDependencies ++ dbmsDependencies ++
    refinedTypeDependenciesList

  val statsDDependencies = List(
    "nl.grons" %% "metrics-scala" % metrics3ScalaV,
    "com.readytalk" % "metrics3-statsd" % metrics3StatsdV
  )

  val gcsFileSystemDependencies = akkaHttpDependencies

  val httpFileSystemDependencies = akkaHttpDependencies

  val ossFileSystemDependencies = googleCloudDependencies ++ aliyunOssDependencies ++ List(
    "com.github.pathikrit" %% "better-files" % betterFilesV
  )

  val statsDProxyDependencies = List(
    "co.fs2" %% "fs2-io" % fs2V,
    "com.iheart" %% "ficus" % ficusV,
    "com.google.cloud" % "google-cloud-nio" % googleCloudNioV
  ) ++ commonDependencies

  val womDependencies = List(
    "com.typesafe.scala-logging" %% "scala-logging" % scalaLoggingV,
    "io.spray" %% "spray-json" % sprayJsonV,
    "org.scalacheck" %% "scalacheck" % scalacheckV % Test,
    "com.github.mpilquist" %% "simulacrum" % simulacrumV,
    "commons-codec" % "commons-codec" % commonsCodecV,
    "eu.timepit" %% "refined" % refinedV
  )

  val wdlDependencies = List(
    "commons-io" % "commons-io" % commonsIoV,
    "org.scala-graph" %% "graph-core" % scalaGraphV,
    "com.chuusai" %% "shapeless" % shapelessV
  ) ++ betterFilesDependencies

  val languageFactoryDependencies = List(
    "com.softwaremill.sttp" %% "core" % sttpV,
    "com.softwaremill.sttp" %% "async-http-client-backend-cats" % sttpV
  )

  val draft2LanguageFactoryDependencies = List(
    "org.mock-server" % "mockserver-netty" % mockserverNettyV % Test
  )

  /*
  The distro artifact contains the actual impl, but transitively includes OSGI bundles that conflict with assembly:
  - https://github.com/owlcs/owlapi/wiki/Documentation/45d8f63d055f820c6ac2ca6c4679a2a7b705449b#howto
  - https://github.com/owlcs/owlapi/issues/455
  - https://github.com/owlcs/owlapi/issues/603

  jcl-over-slf4j.jar is a replacement for commons-logging 1.1.1. Meanwhile our extensive transitive use of apache's
  httpclient has been including commons-logging 1.2 for a while. Now the owl api dependency jcl-over-slf4j is
  conflicting during assembly. As there have been no reported errors AFAIK with commons-logging leaving it in for now.
  However as we use slf4j for cromwell log configuration the correct thing might actually be to exclude commons-logging
  whenever importing httpclient and include jcl-over-slf4j. That way we can control all of our logging in one place.

  - https://www.slf4j.org/legacy.html#jclOverSLF4J
   */
  val owlApiDependencies = List(
    "net.sourceforge.owlapi" % "owlapi-distribution" % owlApiV
      exclude("org.apache.httpcomponents", "httpclient-osgi")
      exclude("org.apache.httpcomponents", "httpcore-osgi")
      exclude("org.slf4j", "jcl-over-slf4j"),
    "org.apache.httpcomponents" % "httpclient-cache" % apacheHttpClientV,
    "org.apache.httpcomponents" % "httpclient" % apacheHttpClientV
  )

  val cwlDependencies = List(
    "com.lihaoyi" %% "ammonite-ops" % ammoniteOpsV,
    "org.broadinstitute" % "heterodon" % heterodonV classifier "single",
    "org.scalactic" %% "scalactic" % scalacticV,
    "org.scalacheck" %% "scalacheck" % scalacheckV % Test,
    "io.circe" %% "circe-optics" % circeOpticsV,
    "org.mozilla" % "rhino" % rhinoV,
    "org.javadelight" % "delight-rhino-sandbox" % delightRhinoSandboxV,
    "org.scalamock" %% "scalamock" % scalamockV % Test,
    "commons-io" % "commons-io" % commonsIoV % Test
  ) ++ circeDependencies ++ womDependencies ++ refinedTypeDependenciesList ++ betterFilesDependencies ++
    owlApiDependencies

  val womtoolDependencies = catsDependencies ++ slf4jBindingDependencies

  val centaurCwlRunnerDependencies = List(
    "com.github.scopt" %% "scopt" % scoptV,
    "io.circe" %% "circe-optics" % circeOpticsV
  ) ++ slf4jBindingDependencies ++ circeDependencies

  val coreDependencies = List(
    "com.google.auth" % "google-auth-library-oauth2-http" % googleOauth2V,
    "com.chuusai" %% "shapeless" % shapelessV,
    "com.github.scopt" %% "scopt" % scoptV,
    "org.scalamock" %% "scalamock" % scalamockV % Test,
  ) ++ akkaStreamDependencies ++ configDependencies ++ catsDependencies ++ circeDependencies ++
    googleApiClientDependencies ++ statsDDependencies ++ betterFilesDependencies ++
    // TODO: We're not using the "F" in slf4j. Core only supports logback, specifically the WorkflowLogger.
    slf4jBindingDependencies

  val databaseMigrationDependencies = liquibaseDependencies ++ dbmsDependencies

  val dockerHashingDependencies = http4sDependencies ++ circeDependencies ++ awsCloudDependencies

  val cromwellApiClientDependencies = List(
    "org.scalaz" %% "scalaz-core" % scalazV,
    "org.typelevel" %% "cats-effect" % catsEffectV,
    "co.fs2" %% "fs2-io" % fs2V % Test,
  ) ++ akkaHttpDependencies ++ betterFilesDependencies ++ catsDependencies

  val centaurDependencies = List(
    "com.github.kxbmap" %% "configs" % configsV,
    "com.google.cloud" % "google-cloud-bigquery" % googleCloudCoreV % IntegrationTest
  ) ++ circeDependencies ++ slf4jBindingDependencies ++ cloudSupportDependencies ++ http4sDependencies

  val engineDependencies = List(
    "commons-codec" % "commons-codec" % commonsCodecV,
    "commons-io" % "commons-io" % commonsIoV,
    "com.storm-enroute" %% "scalameter" % scalameterV
      exclude("com.fasterxml.jackson.core", "jackson-databind")
      exclude("com.fasterxml.jackson.module", "jackson-module-scala")
      exclude("org.scala-tools.testing", "test-interface"),
    "com.fasterxml.jackson.core" % "jackson-databind" % jacksonV,
    "io.github.andrebeat" %% "scala-pool" % scalaPoolV
  ) ++ swaggerUiDependencies ++ akkaHttpDependencies ++ akkaHttpCirceIntegrationDependency ++ circeDependencies

  val serverDependencies = slf4jBindingDependencies

  val cromiamDependencies = List(
    "com.softwaremill.sttp" %% "core" % sttpV,
    "com.softwaremill.sttp" %% "async-http-client-backend-future" % sttpV,
    "com.typesafe.scala-logging" %% "scala-logging" % scalaLoggingV,
    "org.broadinstitute.dsde.workbench" %% "workbench-model" % workbenchModelV,
    "org.broadinstitute.dsde.workbench" %% "workbench-util" % workbenchUtilV
  ) ++ akkaHttpDependencies ++ swaggerUiDependencies ++ slf4jBindingDependencies

  val wes2cromwellDependencies = coreDependencies ++ akkaHttpDependencies

  val backendDependencies = List(
    "org.scalacheck" %% "scalacheck" % scalacheckV % Test,
    "co.fs2" %% "fs2-io" % fs2V
  )

  val bcsBackendDependencies = commonDependencies ++ refinedTypeDependenciesList ++ aliyunBatchComputeDependencies
  val tesBackendDependencies = akkaHttpDependencies
  val sparkBackendDependencies = akkaHttpDependencies

  val testDependencies = List(
    "org.scalatest" %% "scalatest" % scalatestV,
    "org.pegdown" % "pegdown" % pegdownV,
    "org.specs2" %% "specs2-mock" % specs2MockV
  ) ++ slf4jBindingDependencies // During testing, add an slf4j binding for _all_ libraries.

  val kindProjectorPlugin = "org.spire-math" %% "kind-projector" % kindProjectorV
  val paradisePlugin = "org.scalamacros" % "paradise" % paradiseV cross CrossVersion.full

  // Version of the swagger UI to write into config files
  val swaggerUiVersion = swaggerUiV

  val perfDependencies = circeDependencies ++ betterFilesDependencies ++ commonDependencies ++
    googleApiClientDependencies ++ googleCloudDependencies

  val allProjectDependencies =
    backendDependencies ++
      bcsBackendDependencies ++
      centaurCwlRunnerDependencies ++
      centaurDependencies ++
      cloudSupportDependencies ++
      commonDependencies ++
      coreDependencies ++
      cromiamDependencies ++
      cromwellApiClientDependencies ++
      cwlDependencies ++
      databaseMigrationDependencies ++
      databaseSqlDependencies ++
      dockerHashingDependencies ++
      draft2LanguageFactoryDependencies ++
      engineDependencies ++
      gcsFileSystemDependencies ++
      httpFileSystemDependencies ++
      implDrsDependencies ++
      implFtpDependencies ++
      languageFactoryDependencies ++
      ossFileSystemDependencies ++
      perfDependencies ++
      serverDependencies ++
      sparkBackendDependencies ++
      spiDependencies ++
      spiUtilDependencies ++
      statsDProxyDependencies ++
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
   */

  val nettyDependencyOverrides = List(
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

  val rdf4jDependencyOverrides = List(
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

  /*
  If we use a version in one of our projects, that's the one we want all the libraries to use
  ...plus other groups of transitive dependencies shared across multiple projects
   */
  val cromwellDependencyOverrides =
    allProjectDependencies ++
      nettyDependencyOverrides ++
      rdf4jDependencyOverrides
}
