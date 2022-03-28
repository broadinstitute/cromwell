import Publishing._
import Version._

lazy val root = (project in file(".")).
  settings(
    Seq(organization := "org.broadinstitute.cromwell",
    name := "cromwell-client",
    version := createVersion("0.1"),
    scalaVersion := "2.13.8",
    scalacOptions ++= Seq("-feature"),
    compile / javacOptions ++= Seq("-Xlint:deprecation"),
    Compile / packageDoc / publishArtifact := false,
    resolvers += Resolver.mavenLocal,
    updateOptions := updateOptions.value.withGigahorse(false),
    libraryDependencies ++= Seq(
      "io.swagger" % "swagger-annotations" % "1.6.5",
      "com.squareup.okhttp3" % "okhttp" % "4.9.3",
      "com.squareup.okhttp3" % "logging-interceptor" % "4.9.3",
      "com.google.code.gson" % "gson" % "2.9.0",
      "org.apache.commons" % "commons-lang3" % "3.12.0",
      "org.apache.oltu.oauth2" % "org.apache.oltu.oauth2.client" % "1.0.1",
      "javax.ws.rs" % "javax.ws.rs-api" % "2.1.1",
      "javax.annotation" % "javax.annotation-api" % "1.3.2",
      "com.google.code.findbugs" % "jsr305" % "3.0.2",
      "org.threeten" % "threetenbp" % "1.6.0" % Compile,
      "io.gsonfire" % "gson-fire" % "1.8.5" % Compile,
      "junit" % "junit" % "4.13.2" % Test,
      "com.novocode" % "junit-interface" % "0.11" % Test,
      "org.junit.jupiter" % "junit-jupiter-api" % "5.8.2" % Test
    )) ++ publishSettings:_*
  )
