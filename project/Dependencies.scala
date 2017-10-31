import sbt._

object Dependencies {
  val circeV = "0.9.0-M1"

  lazy val akkaV = "2.5.4"
  lazy val akkaHttpV = "10.0.10"

  lazy val slickV = "3.2.0"

  lazy val googleClientApiV = "1.22.0"
  lazy val googleGenomicsServicesApiV = "1.22.0"
  lazy val betterFilesV = "2.17.1"
  lazy val catsV = "1.0.0-MF"
  lazy val mouseV = "0.10-MF"
  lazy val kittensV = "1.0.0-RC0"
  lazy val fs2V = "0.9.7"

  lazy val pegdownV = "1.6.0"
  lazy val scalatestV = "3.0.2"

  // Internal collections of dependencies

  private val fs2Test = "co.fs2" %% "fs2-io" % fs2V % "test"

  private val catsDependencies = List(
    "org.typelevel" %% "cats-core" % catsV,
    "com.github.benhutchison" %% "mouse" % mouseV,
    "org.typelevel" %% "kittens" % kittensV
  ) 

  private val baseDependencies = List(
    "com.iheart" %% "ficus" % "1.4.1",
    "org.scalatest" %% "scalatest" % scalatestV % Test,
    "org.pegdown" % "pegdown" % pegdownV % Test,
    "org.specs2" %% "specs2-mock" % "3.8.9" % Test // 3.9.X doesn't enjoy the spark backend or refined
  ) ++ catsDependencies :+ fs2Test

  private val slf4jBindingDependencies = List(
    // http://logback.qos.ch/dependencies.html
    "ch.qos.logback" % "logback-classic" % "1.2.3",
    "ch.qos.logback" % "logback-access" % "1.2.3",
    "com.getsentry.raven" % "raven-logback" % "8.0.3",
    "org.codehaus.janino" % "janino" % "3.0.7"
  )

  private val slickDependencies = List(
    "com.typesafe.slick" %% "slick" % slickV,
    "com.typesafe.slick" %% "slick-hikaricp" % slickV
  )

  private val liquibaseDependencies = List(
    "org.liquibase" % "liquibase-core" % "3.5.1",
    // This is to stop liquibase from being so noisy by default
    // See: http://stackoverflow.com/questions/20880783/how-to-get-liquibase-to-log-using-slf4j
    "com.mattbertolini" % "liquibase-slf4j" % "2.0.0"
  )

  val akkaHttpDependencies = List(
    "com.typesafe.akka" %% "akka-http" % akkaHttpV,
    "com.typesafe.akka" %% "akka-http-testkit" % akkaHttpV % Test
  )

  val akkaHttpServerDependencies = akkaHttpDependencies :+  "org.webjars" % "swagger-ui" % "3.2.2"

  private val googleApiClientDependencies = List(
    // Used by swagger, but only in tests.  This overrides an older 2.1.3 version of jackson-core brought in by
    // these Google dependencies, but which isn't properly evicted by IntelliJ's sbt integration.
    "com.fasterxml.jackson.core" % "jackson-core" % "2.8.9",
    // The exclusions prevent guava 13 from colliding at assembly time with guava 18 brought in elsewhere.
    "com.google.api-client" % "google-api-client-java6" % googleClientApiV exclude("com.google.guava", "guava-jdk5"),
    "com.google.api-client" % "google-api-client-jackson2" % googleClientApiV exclude("com.google.guava", "guava-jdk5")
  )

  private val googleCloudDependencies = List(
    "io.grpc" % "grpc-core" % "1.5.0",
    "com.google.apis" % "google-api-services-genomics" % ("v1alpha2-rev64-" + googleGenomicsServicesApiV),
    "com.google.cloud" % "google-cloud-nio" % "0.20.1-alpha"
      exclude("com.google.api.grpc", "grpc-google-common-protos")
      exclude("com.google.cloud.datastore", "datastore-v1-protos")
      exclude("org.apache.httpcomponents", "httpclient"),
    "org.apache.httpcomponents" % "httpclient" % "4.5.3"
  )

  private val dbmsDependencies = List(
    "org.hsqldb" % "hsqldb" % "2.3.4",
    /*
    When going to 6.0.x, will need to change the jdbc driver to com.mysql.cj.jdbc.Driver
    - https://dev.mysql.com/doc/connector-j/6.0/en/connector-j-api-changes.html

    The url may also need the parameters:
    - serverTimezone=UTC via http://stackoverflow.com/a/36793896/3320205
    - nullNamePatternMatchesAll=true via https://liquibase.jira.com/browse/CORE-2723
     */
    "mysql" % "mysql-connector-java" % "5.1.42"
  )

  private val refinedTypeDependenciesList = List(
    "eu.timepit" %% "refined" % "0.8.3"
  )

  private val circeYamlDependency = "io.circe" %% "circe-yaml" % "0.7.0-M1"
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

  val gcsFileSystemDependencies = baseDependencies ++ googleApiClientDependencies ++ googleCloudDependencies ++ List (
    "com.github.pathikrit" %% "better-files" % betterFilesV
  )

  val databaseSqlDependencies = baseDependencies ++ slickDependencies ++ dbmsDependencies ++ refinedTypeDependenciesList

  val statsDDependencies = List(
    "nl.grons" %% "metrics-scala" % "3.5.6",
    "com.readytalk" % "metrics3-statsd" % "4.2.0"
  )

  val commonDependencies = List(
    "com.typesafe" % "config" % "1.3.1",
    "org.slf4j" % "slf4j-api" % "1.7.24",
    "com.iheart" %% "ficus" % "1.4.0",
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "org.pegdown" % "pegdown" % pegdownV % Test,
    "com.typesafe.akka" %% "akka-testkit" % akkaV % Test,
    "org.scalatest" %% "scalatest" % scalatestV % Test
  ) ++ catsDependencies ++ slf4jBindingDependencies

  val womDependencies = List(
    "com.typesafe.scala-logging" %% "scala-logging" % "3.6.0",
    "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpV,
    "org.apache.commons" % "commons-lang3" % "3.6",
    "commons-codec" % "commons-codec" % "1.10"
  ) ++ commonDependencies

  val wdlDependencies = List(
    "commons-io" % "commons-io" % "2.5",
    "com.github.pathikrit" %% "better-files" % betterFilesV,
    "org.scala-graph" %% "graph-core" % "1.12.0",
    "com.chuusai" %% "shapeless" % "2.3.2",
    "com.softwaremill.sttp" %% "core" % "0.0.16",
    "com.softwaremill.sttp" %% "async-http-client-backend-cats" % "0.0.16",
    "org.mock-server" % "mockserver-netty" % "3.10.2" % "test"
  ) ++ womDependencies

  val cwlDependencies = List(
    "com.github.pathikrit" %% "better-files" % betterFilesV,
    "com.lihaoyi" %% "ammonite-ops" % "1.0.1",
    "org.typelevel" %% "cats-effect" % "0.4",
    "org.pegdown" % "pegdown" % pegdownV % Test,
    "org.scalactic" %% "scalactic" % "3.0.1",
    "org.scalatest" %% "scalatest" % "3.0.2" % "test",
    "org.scalacheck" %% "scalacheck" % "1.13.4" % "test"
  ) ++ circeDependencies ++ womDependencies ++ refinedTypeDependenciesList

  val womtoolDependencies = wdlDependencies ++ cwlDependencies ++ catsDependencies

  val coreDependencies = List(
    "com.typesafe" % "config" % "1.3.1",
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-slf4j" % akkaV,
    "com.typesafe.akka" %% "akka-testkit" % akkaV % Test,
    "com.google.guava" % "guava" % "22.0",
    "com.google.auth" % "google-auth-library-oauth2-http" % "0.7.0",
    "com.typesafe.akka" %% "akka-stream-testkit" % akkaV,
    "com.chuusai" %% "shapeless" % "2.3.2",
    "com.github.scopt" %% "scopt" % "3.6.0",
    "com.github.pathikrit" %% "better-files" % betterFilesV
  ) ++ baseDependencies ++ googleApiClientDependencies ++ statsDDependencies ++
    // TODO: We're not using the "F" in slf4j. Core only supports logback, specifically the WorkflowLogger.
    slf4jBindingDependencies ++ womDependencies

  val databaseMigrationDependencies = List(
    "com.github.pathikrit" %% "better-files" % betterFilesV % Test
  ) ++ liquibaseDependencies ++ dbmsDependencies

  val cromwellApiClientDependencies = List(
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpV,
    "com.typesafe.akka" %% "akka-stream" % akkaV,
    "com.github.pathikrit" %% "better-files" % betterFilesV,
    "org.scalatest" %% "scalatest" % scalatestV % Test,
    "org.pegdown" % "pegdown" % pegdownV % Test
  ) ++ akkaHttpDependencies

  val engineDependencies = List(
    "commons-codec" % "commons-codec" % "1.10",
    "commons-io" % "commons-io" % "2.5",
    "com.storm-enroute" %% "scalameter" % "0.8.2"
      exclude("com.fasterxml.jackson.core", "jackson-databind")
      exclude("com.fasterxml.jackson.module", "jackson-module-scala")
      exclude("org.scala-tools.testing", "test-interface"),
    "com.fasterxml.jackson.core" % "jackson-databind" % "2.8.9",
    "com.fasterxml.jackson.module" %% "jackson-module-scala" % "2.8.9",
    "io.swagger" % "swagger-parser" % "1.0.22" % Test,
    "org.yaml" % "snakeyaml" % "1.17" % Test
  ) ++ akkaHttpServerDependencies

  val rootDependencies = slf4jBindingDependencies

  val jesBackendDependencies = refinedTypeDependenciesList
  val tesBackendDependencies = akkaHttpDependencies
  val sparkBackendDependencies = akkaHttpDependencies
}
