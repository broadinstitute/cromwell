import com.typesafe.sbt.SbtGit._

name := "CromIam"
organization := "org.broadinstitute"

scalaVersion := "2.12.1"

scalacOptions := List(
  "-Xlint",
  "-feature",
  "-Xmax-classfile-name", "200",
  "-target:jvm-1.8",
  "-encoding", "UTF-8",
  "-unchecked",
  "-deprecation",
  "-Xfuture",
  "-Yno-adapted-args",
  "-Ywarn-dead-code",
  "-Ywarn-numeric-widen",
  "-Ywarn-value-discard",
  "-Ywarn-unused",
  "-Ywarn-unused-import",
  "-Xfatal-warnings"
)

scalacOptions in (Compile, doc) ++= List(
  // http://stackoverflow.com/questions/31488335/scaladoc-2-11-6-fails-on-throws-tag-with-unable-to-find-any-member-to-link#31497874
  "-no-link-warnings"
)

libraryDependencies ++= {
  val akkaV       = "2.4.17"
  val akkaHttpV   = "10.0.6"
  val scalaTestV  = "3.0.1"
  val catsV = "0.9.0"
  val lenthallV = "0.25-903b3c0-SNAP"
  val cromwellV = "27-091ed3b-SNAP"

  val catsDependencies = List(
    "org.typelevel" %% "cats" % catsV
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
  Seq(
    "com.typesafe.akka" %% "akka-actor" % akkaV,
    "com.typesafe.akka" %% "akka-stream" % akkaV,
    "com.typesafe.akka" %% "akka-testkit" % akkaV,
    "com.typesafe.akka" %% "akka-http" % akkaHttpV,
    "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpV,
    "com.typesafe.akka" %% "akka-http-testkit" % akkaHttpV,
    "com.iheart" %% "ficus" % "1.4.0",
    "org.webjars" %  "swagger-ui" % "2.1.1",
    "org.scalatest" %% "scalatest" % scalaTestV % Test,
    "io.swagger" % "swagger-parser" % "1.0.22" % Test,
    "org.yaml" % "snakeyaml" % "1.17" % Test,
    "org.broadinstitute" %% "cromwell-api-client" % cromwellV,
    "org.broadinstitute" %% "lenthall" % lenthallV
  ) ++ catsDependencies
}

imageNames in docker := Seq(
  ImageName(
    namespace = Option("broadinstitute"),
    repository = name.value.toLowerCase,
    tag = git.gitHeadCommit.value)
)

dockerfile in docker := {
  // The assembly task generates a fat JAR file
  val artifact: File = assembly.value
  val artifactTargetPath = s"/app/${artifact.name}"

  new Dockerfile {
    from("openjdk:8")
    expose(8000)
    add(artifact, artifactTargetPath)
    runRaw(s"ln -s $artifactTargetPath /app/cromiam.jar")

    // If you use the 'exec' form for an entry point, shell processing is not performed and
    // environment variable substitution does not occur.  Thus we have to /bin/bash here
    // and pass along any subsequent command line arguments
    // See https://docs.docker.com/engine/reference/builder/#/entrypoint
    entryPoint("/bin/bash", "-c", "java ${JAVA_OPTS} -jar /app/cromiam.jar ${CROMIAM_ARGS} ${*}", "--")
  }
}

buildOptions in docker := BuildOptions(
  cache = false,
  removeIntermediateContainers = BuildOptions.Remove.Always
)

enablePlugins(DockerPlugin)

Revolver.settings
resolvers ++= List(
  "Broad Artifactory Releases" at "https://broadinstitute.jfrog.io/broadinstitute/libs-release/",
  "Broad Artifactory Snapshots" at "https://broadinstitute.jfrog.io/broadinstitute/libs-snapshot/"
)
