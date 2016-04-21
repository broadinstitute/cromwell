import sbt._

object Dependencies {
  lazy val lenthallV = "0.18-e690fc2-SNAPSHOT"
  lazy val wdl4sV = "0.4"
  lazy val sprayV = "1.3.2"
  lazy val DowngradedSprayV = "1.3.1"
  lazy val akkaV = "2.3.12"
  lazy val slickV = "3.1.1"
  lazy val googleClientApiV = "1.20.0"

  val wdl4sDependency = "org.broadinstitute" %% "wdl4s" % wdl4sV

  val testDependencies = List(
    "io.spray" %% "spray-testkit" % sprayV % Test,
    "org.scalatest" %% "scalatest" % "2.2.5" % Test,
    "com.typesafe.akka" %% "akka-testkit" % akkaV % Test,
    "org.yaml" % "snakeyaml" % "1.16" % Test
  )

  val sprayDependencies = List(
    "io.spray" %% "spray-can" % sprayV,
    "io.spray" %% "spray-routing" % sprayV,
    "io.spray" %% "spray-client" % sprayV,
    "io.spray" %% "spray-http" % sprayV,
    "io.spray" %% "spray-json" % DowngradedSprayV
  )

  val googleDependencies = List(
    "com.google.gcloud" % "gcloud-java" % "0.0.9",
    "com.google.api-client" % "google-api-client-java6" % googleClientApiV,
    "com.google.api-client" % "google-api-client-jackson2" % googleClientApiV,
    "com.google.oauth-client" % "google-oauth-client" % googleClientApiV,
    "com.google.cloud.bigdataoss" % "gcsio" % "1.4.4",
    "com.google.apis" % "google-api-services-genomics" % ("v1alpha2-rev14-" + googleClientApiV)
  )

  val gcsFileSystemDependencies = List(
    "org.broadinstitute" %% "lenthall" % lenthallV,
    "org.scalaz" % "scalaz-core_2.11" % "7.1.3"
  ) ++ testDependencies ++ googleDependencies

  val dbDependencies = List(
    "com.typesafe.slick" %% "slick" % slickV,
    "com.typesafe.slick" %% "slick-hikaricp" % slickV,
    "org.hsqldb" % "hsqldb" % "2.3.2",
    "mysql" % "mysql-connector-java" % "5.1.36",
    "org.liquibase" % "liquibase-core" % "3.4.2",
    // This is to stop liquibase from being so noisy by default
    // See: http://stackoverflow.com/questions/20880783/how-to-get-liquibase-to-log-using-slf4j
    "com.mattbertolini" % "liquibase-slf4j" % "2.0.0"
  )

  val coreDependencies = List(
    wdl4sDependency,
    "com.typesafe" % "config" % "1.3.0",
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "org.apache.commons" % "commons-lang3" % "3.4"
  ) ++ testDependencies

  val backendDependencies = coreDependencies

  val engineDependencies = List(
    "org.broadinstitute" %% "lenthall" % lenthallV,
    "com.typesafe.scala-logging" %% "scala-logging" % "3.1.0",
    "org.joda" % "joda-convert" % "1.8.1",
    "org.webjars" % "swagger-ui" % "2.1.1",
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-slf4j" % akkaV,
    "commons-codec" % "commons-codec" % "1.10",
    "commons-io" % "commons-io" % "2.4",
    "ch.qos.logback" % "logback-classic" % "1.1.3",
    "ch.qos.logback" % "logback-access" % "1.1.3",
    "org.codehaus.janino" % "janino" % "2.7.8",
    "org.scalaz" % "scalaz-core_2.11" % "7.1.3",
    "com.github.pathikrit" %% "better-files" % "2.13.0"
  ) ++ coreDependencies ++ sprayDependencies ++ googleDependencies ++ dbDependencies ++ backendDependencies ++ gcsFileSystemDependencies
}
