import sbt._

object Dependencies {
  lazy val lenthallV = "0.25"
  lazy val wdl4sV = "0.14-ac96986-SNAP"

  lazy val akkaV = "2.4.17"
  lazy val akkaHttpV = "10.0.9"

  lazy val slickV = "3.2.0"

  lazy val googleClientApiV = "1.22.0"
  lazy val googleGenomicsServicesApiV = "1.22.0"
  lazy val betterFilesV = "2.17.1"
  lazy val catsV = "0.9.0"
  lazy val fs2V = "0.9.7"

  lazy val pegdownV = "1.6.0"
  lazy val scalatestV = "3.0.2"

  // Internal collections of dependencies

  private val fs2Test = "co.fs2" %% "fs2-io" % fs2V % "test"

  private val catsDependencies = List(
    "org.typelevel" %% "cats" % catsV,
    "com.github.benhutchison" %% "mouse" % "0.9"
  ) map (_
    /*
    Exclude test framework cats-laws and its transitive dependency scalacheck.
    If sbt detects scalacheck, it tries to run it.
    Explicitly excluding the two problematic artifacts instead of including the three (or four?).
    https://github.com/typelevel/cats/tree/v0.7.2#getting-started
    Re "_2.12", see also: https://github.com/sbt/sbt/issues/1518
     */
    exclude("org.typelevel", "cats-laws_2.12")
    exclude("org.typelevel", "cats-kernel-laws_2.12")
    )

  private val baseDependencies = List(
    "org.broadinstitute" %% "lenthall" % lenthallV,
    "com.iheart" %% "ficus" % "1.4.1",
    "org.scalatest" %% "scalatest" % scalatestV % Test,
    "org.pegdown" % "pegdown" % pegdownV % Test,
    "org.specs2" %% "specs2-mock" % "3.8.9" % Test // 3.9.X doesn't enjoy the spark backend or refined
  ) ++ catsDependencies :+ fs2Test

  private val slf4jBindingDependencies = List(
    // http://logback.qos.ch/dependencies.html
    "ch.qos.logback" % "logback-classic" % "1.2.3",
    "ch.qos.logback" % "logback-access" % "1.2.3",
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

  val akkaHttpServerDependencies = akkaHttpDependencies :+  "org.webjars" % "swagger-ui" % "2.1.1"

  private val googleApiClientDependencies = List(
    // Used by swagger, but only in tests.  This overrides an older 2.1.3 version of jackson-core brought in by
    // these Google dependencies, but which isn't properly evicted by IntelliJ's sbt integration.
    "com.fasterxml.jackson.core" % "jackson-core" % "2.8.9",
    // The exclusions prevent guava 13 from colliding at assembly time with guava 18 brought in elsewhere.
    "com.google.api-client" % "google-api-client-java6" % googleClientApiV exclude("com.google.guava", "guava-jdk5"),
    "com.google.api-client" % "google-api-client-jackson2" % googleClientApiV exclude("com.google.guava", "guava-jdk5")
  )

  private val googleCloudDependencies = List(
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
    "org.scala-lang" % "scala-compiler" % Settings.ScalaVersion,
    "eu.timepit" %% "refined" % "0.8.2"
  )

  // Sub-project dependencies, added in addition to any dependencies inherited from .dependsOn().

  val gcsFileSystemDependencies = baseDependencies ++ googleApiClientDependencies ++ googleCloudDependencies ++ List (
    "com.github.pathikrit" %% "better-files" % betterFilesV
  )

  val databaseSqlDependencies = baseDependencies ++ slickDependencies ++ dbmsDependencies ++ refinedTypeDependenciesList

  val coreDependencies = List(
    "com.typesafe.scala-logging" %% "scala-logging" % "3.6.0",
    "org.broadinstitute" %% "wdl4s-wdl" % wdl4sV,
    "org.apache.commons" % "commons-lang3" % "3.6",
    "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpV,
    "com.typesafe" % "config" % "1.3.1",
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-slf4j" % akkaV,
    "com.typesafe.akka" %% "akka-testkit" % akkaV % Test,
    "com.google.guava" % "guava" % "22.0",
    "com.google.auth" % "google-auth-library-oauth2-http" % "0.7.0",
    "com.typesafe.akka" %% "akka-stream-testkit" % akkaV,
    "com.chuusai" %% "shapeless" % "2.3.2",
    "com.github.scopt" %% "scopt" % "3.6.0"
  ) ++ baseDependencies ++ googleApiClientDependencies ++
    // TODO: We're not using the "F" in slf4j. Core only supports logback, specifically the WorkflowLogger.
    slf4jBindingDependencies

  val databaseMigrationDependencies = List(
    "com.github.pathikrit" %% "better-files" % betterFilesV % Test
  ) ++ liquibaseDependencies ++ dbmsDependencies

  val cromwellApiClientDependencies = List(
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-http" % akkaHttpV,
    "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpV,
    "com.github.pathikrit" %% "better-files" % betterFilesV,
    "org.scalatest" %% "scalatest" % scalatestV % Test,
    "org.pegdown" % "pegdown" % pegdownV % Test
  )

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
