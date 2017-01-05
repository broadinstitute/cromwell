import sbt._

object Dependencies {
  lazy val lenthallV = "0.20"
  lazy val wdl4sV = "0.8"
  lazy val sprayV = "1.3.3"
  /*
  spray-json is an independent project from the "spray suite"
  - https://github.com/spray/spray
  - https://github.com/spray/spray-json
  - http://spray.io/documentation/1.2.2/spray-httpx/spray-json-support/
  - http://doc.akka.io/docs/akka/2.4/scala/http/common/json-support.html#akka-http-spray-json
   */
  lazy val sprayJsonV = "1.3.2"
  lazy val akkaV = "2.4.14"
  lazy val slickV = "3.1.1"
  lazy val googleClientApiV = "1.22.0"
  lazy val googleGenomicsServicesApiV = "1.20.0"
  lazy val betterFilesV = "2.16.0"
  lazy val catsV = "0.7.2"

  // Internal collections of dependencies

  private val catsDependencies = List(
    "org.typelevel" %% "cats" % "0.7.2",
    "com.github.benhutchison" %% "mouse" % "0.5"
  ) map (_
    /*
    Exclude test framework cats-laws and its transitive dependency scalacheck.
    If sbt detects scalacheck, it tries to run it.
    Explicitly excluding the two problematic artifacts instead of including the three (or four?).
    https://github.com/typelevel/cats/tree/v0.7.2#getting-started
    Re "_2.11", see also: https://github.com/sbt/sbt/issues/1518
     */
    exclude("org.typelevel", "cats-laws_2.11")
    exclude("org.typelevel", "cats-kernel-laws_2.11")
    )

  private val baseDependencies = List(
    "org.broadinstitute" %% "lenthall" % lenthallV,
    "com.iheart" %% "ficus" % "1.3.0",
    "org.scalatest" %% "scalatest" % "3.0.0" % Test,
    "org.pegdown" % "pegdown" % "1.6.0" % Test,
    "org.specs2" %% "specs2-mock" % "3.8.5" % Test
  ) ++ catsDependencies

  private val slf4jBindingDependencies = List(
    // http://logback.qos.ch/dependencies.html
    "ch.qos.logback" % "logback-classic" % "1.1.7",
    "ch.qos.logback" % "logback-access" % "1.1.7",
    "org.codehaus.janino" % "janino" % "3.0.1"
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

  private val sprayServerDependencies = List(
    "io.spray" %% "spray-can" % sprayV,
    "io.spray" %% "spray-routing-shapeless2" % sprayV,
    "io.spray" %% "spray-http" % sprayV,
    "io.spray" %% "spray-testkit" % sprayV % Test
  )

  private val googleApiClientDependencies = List(
    // Used by swagger, but only in tests.  This overrides an older 2.1.3 version of jackson-core brought in by
    // these Google dependencies, but which isn't properly evicted by IntelliJ's sbt integration.
    "com.fasterxml.jackson.core" % "jackson-core" % "2.8.2",
    // The exclusions prevent guava 13 from colliding at assembly time with guava 18 brought in elsewhere.
    "com.google.api-client" % "google-api-client-java6" % googleClientApiV exclude("com.google.guava", "guava-jdk5"),
    "com.google.api-client" % "google-api-client-jackson2" % googleClientApiV exclude("com.google.guava", "guava-jdk5")
  )

  private val googleCloudDependencies = List(
    "com.google.apis" % "google-api-services-genomics" % ("v1alpha2-rev14-" + googleGenomicsServicesApiV),
    "com.google.cloud" % "google-cloud-nio" % "0.3.0"
      exclude("com.google.api.grpc", "grpc-google-common-protos")
      exclude("com.google.cloud.datastore", "datastore-v1-protos")
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
    "mysql" % "mysql-connector-java" % "5.1.39"
  )

  // Sub-project dependencies, added in addition to any dependencies inherited from .dependsOn().

  val gcsFileSystemDependencies = baseDependencies ++ googleApiClientDependencies ++ googleCloudDependencies ++ List (
    "com.github.pathikrit" %% "better-files" % betterFilesV
  )

  val databaseSqlDependencies = baseDependencies ++ slickDependencies ++ dbmsDependencies

  val coreDependencies = List(
    "com.typesafe.scala-logging" %% "scala-logging" % "3.4.0",
    "org.broadinstitute" %% "wdl4s" % wdl4sV,
    "org.apache.commons" % "commons-lang3" % "3.4",
    "io.spray" %% "spray-json" % sprayJsonV,
    "com.typesafe" % "config" % "1.3.0",
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-slf4j" % akkaV,
    "com.typesafe.akka" %% "akka-testkit" % akkaV % Test,
    "com.google.guava" % "guava" % "20.0"
  ) ++ baseDependencies ++ googleApiClientDependencies ++
    // TODO: We're not using the "F" in slf4j. Core only supports logback, specifically the WorkflowLogger.
    slf4jBindingDependencies

  val databaseMigrationDependencies = List(
    "com.github.pathikrit" %% "better-files" % betterFilesV % Test
  ) ++ liquibaseDependencies ++ dbmsDependencies

  val htCondorBackendDependencies = List(
    "com.twitter" %% "chill" % "0.8.0",
    "org.mongodb" %% "casbah" % "3.0.0"
  )

  val sparkBackendDependencies = List(
    "io.spray" %% "spray-client" % sprayV
  ) ++ sprayServerDependencies

  val engineDependencies = List(
    "org.webjars" % "swagger-ui" % "2.1.1",
    "commons-codec" % "commons-codec" % "1.10",
    "commons-io" % "commons-io" % "2.5",
    "io.swagger" % "swagger-parser" % "1.0.22" % Test,
    "org.yaml" % "snakeyaml" % "1.17" % Test
  ) ++ sprayServerDependencies

  val rootDependencies = slf4jBindingDependencies
}
