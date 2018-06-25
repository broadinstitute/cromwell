import sbt._

object Dependencies {
  private val akkaHttpV = "10.0.10"
  private val akkaV = "2.5.4"
  private val alibabaCloudCoreV = "3.6.0"
  private val alibabaCloudOssV = "3.1.0"
  private val alibabaCloudBcsV = "5.3.2"
  private val ammoniteOpsV = "1.0.1"
  private val apacheHttpClientV = "4.5.3"
  private val apacheHttpCoreV = "4.4.6"
  private val awsSdkV = "2.0.0-preview-9"
  private val s3fsV = "1.0.1"
  private val betterFilesV = "2.17.1"
  private val catsEffectV = "0.10"
  private val catsV = "1.0.1"
  private val circeV = "0.9.3"
  private val circeYamlV = "0.7.0"
  private val commonsCodecV = "1.10"
  private val commonsIoV = "2.5"
  private val commonsLang3V = "3.5"
  private val commonsLoggingV = "1.2"
  private val commonsTextV = "1.1"
  private val configsV = "0.4.4"
  private val delightRhinoSandboxV = "0.0.8"
  private val errorProneAnnotationsV = "2.0.19"
  private val ficusV = "1.4.1"
  private val fs2V = "0.10.2"
  private val gaxV = "1.9.0"
  private val googleApiClientV = "1.23.0"
  private val googleCloudComputeV = "0.26.0-alpha"
  private val googleCloudCoreV = "1.8.0"
  private val googleCloudNioV = "0.20.1-alpha"
  private val googleCredentialsV = "0.8.0"
  private val googleGenomicsServicesV2ApiV = "v2alpha1-rev15-1.23.0"
  private val googleGenomicsServicesV1ApiV = "v1alpha2-rev495-1.23.0"
  private val googleHttpClientV = googleApiClientV
  private val googleOauth2V = "0.8.0"
  private val googleOauthClientV = googleApiClientV
  private val grpcV = "1.5.0"
  private val guavaV = "22.0"
  private val hsqldbV = "2.3.4"
  private val jacksonV = "2.9.4"
  private val janinoV = "3.0.7"
  private val jodaTimeV = "2.9.4"
  private val jsr305V = "3.0.0"
  private val kindProjectorV = "0.9.4"
  private val kittensV = "1.0.0-RC3"
  private val liquibaseSlf4jV = "2.0.0"
  private val liquibaseV = "3.5.1"
  private val logbackV = "1.2.3"
  private val metrics3StatsdV = "4.2.0"
  private val metricsScalaV = "3.5.6"
  private val mockserverNettyV = "3.10.2"
  private val mouseV = "0.10-MF"
  private val mysqlV = "5.1.42"
  private val nettyHandlerV = "4.1.22.Final"
  private val owlApiV = "5.1.4"
  private val paradiseV = "2.1.0"
  private val pegdownV = "1.6.0"
  private val protoGoogleCommonProtosV = "0.1.21"
  private val protoGoogleIamV1V = "0.1.21"
  private val protobufJavaV = "3.3.1"
  private val ravenLogbackV = "8.0.3"
  private val reactiveStreamsV = "1.0.1"
  private val refinedV = "0.8.3"
  private val rhinoV = "1.7.8"
  private val scalaGraphV = "1.12.0"
  private val scalaLoggingV = "3.7.1"
  private val scalaPoolV = "0.4.1"
  private val scalaXmlV = "1.0.6"
  private val scalacheckV = "1.13.4"
  private val scalacticV = "3.0.1"
  private val scalameterV = "0.8.2"
  private val scalamockV = "4.0.0"
  private val scalatestV = "3.0.2"
  private val scalazV = "7.2.17"
  private val scoptV = "3.6.0"
  private val shapelessV = "2.3.3"
  private val simulacrumV = "0.12.0"
  private val slf4jV = "1.7.24"
  private val slickV = "3.2.3"
  private val slickCatsV = "0.7.1"
  private val snakeyamlV = "1.17"
  private val specs2MockV = "3.8.9" // 3.9.X doesn't enjoy the spark backend or refined
  private val sprayJsonV = "1.3.3"
  private val squantV = "1.3.0"
  private val sttpV = "1.1.12"
  private val swaggerParserV = "1.0.22"
  private val swaggerUiV = "3.2.2"
  private val typesafeConfigV = "1.3.1"
  private val workbenchGoogleV = "0.15-2fc79a3"
  private val workbenchModelV = "0.10-6800f3a"
  private val workbenchUtilV = "0.3-f3ce961"

  /*
  If you see warnings from SBT about evictions, insert a specific dependency version into this list.
   */
  val cromwellDependencyOverrides = List(
    "ch.qos.logback" % "logback-classic" % logbackV,
    "ch.qos.logback" % "logback-core" % logbackV,
    "com.fasterxml.jackson.core" % "jackson-annotations" % jacksonV,
    "com.fasterxml.jackson.core" % "jackson-core" % jacksonV,
    "com.fasterxml.jackson.module" %% "jackson-module-scala" % jacksonV,
    "com.google.api" % "gax" % gaxV,
    "com.google.api-client" % "google-api-client" % googleApiClientV,
    "com.google.api.grpc" % "proto-google-common-protos" % protoGoogleCommonProtosV,
    "com.google.api.grpc" % "proto-google-iam-v1" % protoGoogleIamV1V,
    "com.google.auth" % "google-auth-library-credentials" % googleCredentialsV,
    "com.google.auth" % "google-auth-library-oauth2-http" % googleOauth2V,
    "com.google.cloud" % "google-cloud-core" % googleCloudCoreV,
    "com.google.cloud" % "google-cloud-core-http" % googleCloudCoreV,
    "com.google.code.findbugs" % "jsr305" % jsr305V,
    "com.google.errorprone" % "error_prone_annotations" % errorProneAnnotationsV,
    "com.google.guava" % "guava" % guavaV,
    "com.google.http-client" % "google-http-client" % googleHttpClientV,
    "com.google.http-client" % "google-http-client-appengine" % googleHttpClientV,
    "com.google.http-client" % "google-http-client-jackson" % googleHttpClientV,
    "com.google.http-client" % "google-http-client-jackson2" % googleHttpClientV,
    "com.google.oauth-client" % "google-oauth-client" % googleOauthClientV,
    "com.google.protobuf" % "protobuf-java" % protobufJavaV,
    "com.google.protobuf" % "protobuf-java-util" % protobufJavaV,
    "com.typesafe" % "config" % typesafeConfigV,
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-stream" % akkaV,
    "commons-codec" % "commons-codec" % commonsCodecV,
    "commons-io" % "commons-io" % commonsIoV,
    "commons-logging" % "commons-logging" % commonsLoggingV,
    "eu.timepit" %% "refined" % refinedV,
    "io.grpc" % "grpc-context" % grpcV,
    "io.netty" % "netty-handler" % nettyHandlerV,
    "io.spray" %% "spray-json" % sprayJsonV,
    "joda-time" % "joda-time" % jodaTimeV,
    "org.apache.commons" % "commons-lang3" % commonsLang3V,
    "org.apache.httpcomponents" % "httpclient" % apacheHttpClientV,
    "org.apache.httpcomponents" % "httpcore" % apacheHttpCoreV,
    "org.reactivestreams" % "reactive-streams" % reactiveStreamsV,
    "org.scala-lang.modules" %% "scala-xml" % scalaXmlV,
    "org.slf4j" % "slf4j-api" % slf4jV,
    "org.typelevel" %% "cats-core" % catsV,
    "org.typelevel" %% "cats-kernel" % catsV,
    "org.yaml" % "snakeyaml" % snakeyamlV
  )

  // Internal collections of dependencies

  private val betterFilesDependencies = List(
    "com.github.pathikrit" %% "better-files" % betterFilesV
  )

  private val catsDependencies = List(
    "org.typelevel" %% "cats-core" % catsV,
    "org.typelevel" %% "alleycats-core" % catsV,
    "com.github.benhutchison" %% "mouse" % mouseV,
    "org.typelevel" %% "kittens" % kittensV
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
    "ch.qos.logback" % "logback-classic" % logbackV,
    "ch.qos.logback" % "logback-access" % logbackV,
    "com.typesafe.akka" %% "akka-slf4j" % akkaV,
    "com.getsentry.raven" % "raven-logback" % ravenLogbackV,
    "org.codehaus.janino" % "janino" % janinoV
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

  private val akkaHttpDependencies = List(
    "com.typesafe.akka" %% "akka-http" % akkaHttpV,
    "com.typesafe.akka" %% "akka-http-testkit" % akkaHttpV % Test
  )

  private val swaggerUiDependencies = List(
    "org.webjars" % "swagger-ui" % swaggerUiV,
    "io.swagger" % "swagger-parser" % swaggerParserV % Test,
    "org.yaml" % "snakeyaml" % snakeyamlV % Test
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

  // The v1 dependency has been cloned in the broad artifactory so that we can have the 2 versions co-exist in the same jar
  private val googleGenomicsV1Dependency = List(
    "org.broadinstitute" % "cromwell-google-api-services-genomics" % googleGenomicsServicesV1ApiV
      exclude("com.google.guava", "guava-jdk5")
  )

  private val googleGenomicsV2Dependency = List(
    "com.google.apis" % "google-api-services-genomics" % googleGenomicsServicesV2ApiV
      exclude("com.google.guava", "guava-jdk5")
  )

  private val awsCloudDependencies = List(
    "com.fasterxml.jackson.core" % "jackson-annotations" % jacksonV,
    "software.amazon.awssdk" % "aws-sdk-java" % awsSdkV,
    "org.lerch" % "s3fs" % s3fsV
      exclude("org.slf4j", "jcl-over-slf4j")
  )

  private val googleCloudDependencies = List(
    "io.grpc" % "grpc-core" % grpcV,
    "com.google.guava" % "guava" % guavaV,
    "com.google.cloud" % "google-cloud-nio" % googleCloudNioV
      exclude("com.google.api.grpc", "grpc-google-common-protos")
      exclude("com.google.cloud.datastore", "datastore-v1-protos")
      exclude("org.apache.httpcomponents", "httpclient"),
    "com.google.cloud" % "google-cloud-compute" % googleCloudComputeV,
    "org.broadinstitute.dsde.workbench" %% "workbench-google" % workbenchGoogleV
      exclude("com.google.apis", "google-api-services-genomics"),
    "org.apache.httpcomponents" % "httpclient" % apacheHttpClientV
  ) ++ googleGenomicsV1Dependency ++ googleGenomicsV2Dependency

  private val aliyunOssDependencies = List(
    "com.aliyun.oss" % "aliyun-sdk-oss" % alibabaCloudOssV
      // stax is included twice by oss 3.1.0 and cause assembly merge conflicts via stax vs. javax.xml.stream
      exclude("stax", "stax-api")
  )

  private val aliyunBatchComputeDependencies = List(
    "com.aliyun" % "aliyun-java-sdk-core" % alibabaCloudCoreV,
    "com.aliyun" % "aliyun-java-sdk-batchcompute" % alibabaCloudBcsV
  )

  private val dbmsDependencies = List(
    "org.hsqldb" % "hsqldb" % hsqldbV,
    /*
    When going to 6.0.x, will need to change the jdbc driver to com.mysql.cj.jdbc.Driver
    - https://dev.mysql.com/doc/connector-j/6.0/en/connector-j-api-changes.html

    The url may also need the parameters:
    - serverTimezone=UTC via http://stackoverflow.com/a/36793896/3320205
    - nullNamePatternMatchesAll=true via https://liquibase.jira.com/browse/CORE-2723
     */
    "mysql" % "mysql-connector-java" % mysqlV
  )

  private val refinedTypeDependenciesList = List(
    "eu.timepit" %% "refined" % refinedV
  )

  private val circeYamlDependency = "io.circe" %% "circe-yaml" % circeYamlV

  private val circeDependencies = List(
    "core",
    "parser",
    "generic",
    "generic-extras",
    "shapes",
    "refined",
    "literal"
  ).map(m => "io.circe" %% s"circe-$m" % circeV) :+ circeYamlDependency

  // Sub-project dependencies, added in addition to any dependencies inherited from .dependsOn().

  val cloudSupportDependencies = googleApiClientDependencies ++ googleCloudDependencies ++ betterFilesDependencies ++ awsCloudDependencies

  val databaseSqlDependencies = configDependencies ++ catsDependencies ++ slickDependencies ++ dbmsDependencies ++
    refinedTypeDependenciesList

  val statsDDependencies = List(
    "nl.grons" %% "metrics-scala" % metricsScalaV,
    "com.readytalk" % "metrics3-statsd" % metrics3StatsdV
  )

  val ossFileSystemDependencies = googleCloudDependencies ++ aliyunOssDependencies ++ List (
    "com.github.pathikrit" %% "better-files" % betterFilesV
  )

  val commonDependencies = List(
    "org.slf4j" % "slf4j-api" % slf4jV,
    "org.typelevel" %% "cats-effect" % catsEffectV,
    "org.apache.commons" % "commons-lang3" % commonsLang3V
  ) ++ catsDependencies ++ configDependencies

  val womDependencies = List(
    "com.typesafe.scala-logging" %% "scala-logging" % scalaLoggingV,
    "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpV,
    "org.scalacheck" %% "scalacheck" % scalacheckV % Test,
    "org.apache.commons" % "commons-text" % commonsTextV,
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
      exclude("org.apache.httpcomponents", "httpcore-osgi")
      exclude("org.slf4j", "jcl-over-slf4j"),
    "org.apache.httpcomponents" % "httpclient" % apacheHttpClientV
  )

  val cwlDependencies = List(
    "com.lihaoyi" %% "ammonite-ops" % ammoniteOpsV,
    "org.scalactic" %% "scalactic" % scalacticV,
    "org.scalacheck" %% "scalacheck" % scalacheckV % Test,
    "io.circe" %% "circe-optics" % circeV,
    "org.mozilla" % "rhino" % rhinoV,
    "org.javadelight" % "delight-rhino-sandbox" % delightRhinoSandboxV,
    "org.scalamock" %% "scalamock" % scalamockV % Test,
    "commons-io" % "commons-io" % commonsIoV % Test
  ) ++ circeDependencies ++ womDependencies ++ refinedTypeDependenciesList ++ betterFilesDependencies ++
    owlApiDependencies

  val womtoolDependencies = catsDependencies ++ slf4jBindingDependencies

  val centaurCwlRunnerDependencies = List(
    "com.github.scopt" %% "scopt" % scoptV,
    "io.circe" %% "circe-optics" % circeV
  ) ++ slf4jBindingDependencies ++ circeDependencies

  val coreDependencies = List(
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-testkit" % akkaV % Test,
    "com.typesafe.akka" %% "akka-stream" % akkaV,
    "com.typesafe.akka" %% "akka-stream-testkit" % akkaV % Test,
    "com.google.auth" % "google-auth-library-oauth2-http" % googleOauth2V,
    "com.chuusai" %% "shapeless" % shapelessV,
    "com.github.scopt" %% "scopt" % scoptV,
    "org.typelevel"  %% "squants"  % squantV
  ) ++ configDependencies ++ catsDependencies ++ googleApiClientDependencies ++ statsDDependencies ++
    betterFilesDependencies ++
    // TODO: We're not using the "F" in slf4j. Core only supports logback, specifically the WorkflowLogger.
    slf4jBindingDependencies

  val databaseMigrationDependencies = liquibaseDependencies ++ dbmsDependencies

  val cromwellApiClientDependencies = List(
    "org.scalaz" %% "scalaz-core" % scalazV,
    "co.fs2" %% "fs2-io" % fs2V % Test,
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpV,
    "com.typesafe.akka" %% "akka-stream" % akkaV
  ) ++ akkaHttpDependencies ++ betterFilesDependencies

  val centaurDependencies = List(
    "com.github.kxbmap" %% "configs" % configsV
  ) ++ circeDependencies ++ slf4jBindingDependencies ++ cloudSupportDependencies

  val engineDependencies = List(
    "commons-codec" % "commons-codec" % commonsCodecV,
    "commons-io" % "commons-io" % commonsIoV,
    "com.storm-enroute" %% "scalameter" % scalameterV
      exclude("com.fasterxml.jackson.core", "jackson-databind")
      exclude("com.fasterxml.jackson.module", "jackson-module-scala")
      exclude("org.scala-tools.testing", "test-interface"),
    "com.fasterxml.jackson.core" % "jackson-databind" % jacksonV,
    "io.github.andrebeat" %% "scala-pool" % scalaPoolV
  ) ++ swaggerUiDependencies ++ akkaHttpDependencies ++ circeDependencies

  val serverDependencies = slf4jBindingDependencies

  val cromiamDependencies = List(
    "com.softwaremill.sttp" %% "core" % sttpV,
    "com.softwaremill.sttp" %% "async-http-client-backend-future" % sttpV,
    "com.typesafe.scala-logging" %% "scala-logging" % scalaLoggingV,
    "org.broadinstitute.dsde.workbench" %% "workbench-model" % workbenchModelV,
    "org.broadinstitute.dsde.workbench" %% "workbench-util" % workbenchUtilV
  ) ++ akkaHttpDependencies ++ catsDependencies ++ swaggerUiDependencies

  val backendDependencies = List(
    "org.scalacheck" %% "scalacheck" % scalacheckV % Test,
    "co.fs2" %% "fs2-io" % fs2V % Test
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
}
