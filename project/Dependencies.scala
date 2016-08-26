import sbt._

object Dependencies {
  lazy val lenthallV = "0.18-e690fc2-SNAPSHOT"
  lazy val wdl4sV = "0.5-ed015a6-SNAPSHOT"
  lazy val sprayV = "1.3.2"
  /*
  spray-json is an independent project from the "spray suite"
  - https://github.com/spray/spray
  - https://github.com/spray/spray-json
  - http://spray.io/documentation/1.2.2/spray-httpx/spray-json-support/
  - http://doc.akka.io/docs/akka/2.4/scala/http/common/json-support.html#akka-http-spray-json
   */
  lazy val sprayJsonV = "1.3.1"
  lazy val akkaV = "2.3.15"
  lazy val slickV = "3.1.1"
  lazy val googleClientApiV = "1.20.0"
  lazy val betterFilesV = "2.13.0"
  lazy val scalazCoreV = "7.1.3"

  // Internal collections of dependencies

  private val baseDependencies = List(
    "org.broadinstitute" %% "lenthall" % lenthallV,
    "org.scalaz" %% "scalaz-core" % scalazCoreV,
    "org.scalatest" %% "scalatest" % "3.0.0" % Test,
    "org.specs2" %% "specs2" % "2.3.13" % Test
  )

  private val slf4jBindingDependencies: List[ModuleID] = List(
    // http://logback.qos.ch/dependencies.html
    "ch.qos.logback" % "logback-classic" % "1.1.3",
    "ch.qos.logback" % "logback-access" % "1.1.3",
    "org.codehaus.janino" % "janino" % "2.7.8"
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
    "com.fasterxml.jackson.core" % "jackson-core" % "2.4.5",
    // The exclusions prevent guava 13 from colliding at assembly time with guava 18 brought in elsewhere.
    "com.google.api-client" % "google-api-client-java6" % googleClientApiV exclude("com.google.guava", "guava-jdk5"),
    "com.google.api-client" % "google-api-client-jackson2" % googleClientApiV exclude("com.google.guava", "guava-jdk5")
  )

  private val googleCloudDependencies = List(
    "com.google.gcloud" % "gcloud-java" % "0.0.9",
    "com.google.oauth-client" % "google-oauth-client" % googleClientApiV,
    "com.google.cloud.bigdataoss" % "gcsio" % "1.4.4",
    "com.google.apis" % "google-api-services-genomics" % ("v1alpha2-rev14-" + googleClientApiV)
  )

  // Sub-project dependencies, added in addition to any dependencies inherited from .dependsOn().

  val gcsFileSystemDependencies = baseDependencies ++ googleApiClientDependencies ++ googleCloudDependencies

  val databaseDependencies = List(
    "com.typesafe.slick" %% "slick" % slickV,
    "com.typesafe.slick" %% "slick-hikaricp" % slickV,
    "org.hsqldb" % "hsqldb" % "2.3.2",
    "mysql" % "mysql-connector-java" % "5.1.36",
    "com.github.pathikrit" %% "better-files" % betterFilesV % Test
  ) ++ baseDependencies ++ liquibaseDependencies

  val coreDependencies = List(
    "com.typesafe.scala-logging" %% "scala-logging" % "3.1.0",
    "org.broadinstitute" %% "wdl4s" % wdl4sV,
    "org.apache.commons" % "commons-lang3" % "3.4",
    "io.spray" %% "spray-json" % sprayJsonV,
    "com.typesafe" % "config" % "1.3.0",
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-slf4j" % akkaV,
    "com.typesafe.akka" %% "akka-testkit" % akkaV % Test
  ) ++ baseDependencies ++ googleApiClientDependencies ++
    // TODO: We're not using the "F" in slf4j. Core only supports logback, specifically the WorkflowLogger.
    slf4jBindingDependencies

  val htCondorBackendDependencies = List(
    "com.twitter" %% "chill" % "0.8.0",
    "org.mongodb" %% "casbah" % "3.0.0"
  )

  val engineDependencies = List(
    "com.typesafe.scala-logging" %% "scala-logging" % "3.1.0",
    "org.webjars" % "swagger-ui" % "2.1.1",
    "commons-codec" % "commons-codec" % "1.10",
    "commons-io" % "commons-io" % "2.4",
    "org.scalaz" %% "scalaz-core" % scalazCoreV,
    "com.github.pathikrit" %% "better-files" % betterFilesV,
    "io.swagger" % "swagger-parser" % "1.0.19" % Test,
    "org.yaml" % "snakeyaml" % "1.16" % Test
  ) ++ sprayServerDependencies

  val rootDependencies = slf4jBindingDependencies
}
