import java.io.FileNotFoundException

import Version.cromwellVersion
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtdocker.DockerPlugin.autoImport._
import ContinuousIntegration._

import scala.sys.process._

object Publishing {

  val dockerTags = settingKey[Seq[String]]("The tags for docker builds.")
  val dockerPushCheck = taskKey[Unit]("Check that the public repository exists in DockerHub.")

  val dockerSettings: Seq[Setting[_]] = List(
    /*
    NOTE: Gave up fighting with SBT settings. Using an environment variable instead.

    The below "just works", assuming womtool docker image building is also enabled, setting the right image names and
    versions.

    `sbt 'show docker::imageNames'` returns:
      ArrayBuffer(broadinstitute/womtool:30-c33be41-SNAP)
      ArrayBuffer(broadinstitute/cromwell:30-c33be41-SNAP)

    `CROMWELL_SBT_DOCKER_TAGS=dev,develop sbt 'show docker::imageNames'` returns:
      ArrayBuffer(broadinstitute/womtool:dev, broadinstitute/womtool:develop)
      ArrayBuffer(broadinstitute/cromwell:dev, broadinstitute/cromwell:develop)
    */
    dockerTags := {
      val versionsCsv = if (Version.isSnapshot) version.value else s"$cromwellVersion,${version.value}"
      sys.env.getOrElse("CROMWELL_SBT_DOCKER_TAGS", versionsCsv).split(",")
    },
    imageNames in docker := dockerTags.value map { tag =>
      ImageName(namespace = Option("broadinstitute"), repository = name.value, tag = Option(tag))
    },
    dockerfile in docker := {
      // The assembly task generates a fat JAR file
      val artifact: File = assembly.value
      val artifactTargetPath = s"/app/${artifact.name}"
      val projectName = name.value

      new Dockerfile {
        from("openjdk:8")
        expose(8000)
        add(artifact, artifactTargetPath)
        runRaw(s"ln -s $artifactTargetPath /app/$projectName.jar")

        // If you use the 'exec' form for an entry point, shell processing is not performed and
        // environment variable substitution does not occur.  Thus we have to /bin/bash here
        // and pass along any subsequent command line arguments
        // See https://docs.docker.com/engine/reference/builder/#/entrypoint
        entryPoint(
          "/bin/bash",
          "-c",
          s"java $${JAVA_OPTS} -jar /app/$projectName.jar $${${projectName.toUpperCase.replaceAll("-", "_")}_ARGS} $${*}",
          "--"
        )
      }
    },
    buildOptions in docker := BuildOptions(
      cache = false,
      removeIntermediateContainers = BuildOptions.Remove.Always
    ),
  )

  def dockerPushSettings(pushEnabled: Boolean): Seq[Setting[_]] = {
    if (pushEnabled) {
      List(
        dockerPushCheck := {
          val projectName = name.value
          val repositoryName = s"broadinstitute/$projectName"
          val repositoryUrl = s"https://registry.hub.docker.com/v2/repositories/$repositoryName/"
          try {
            url(repositoryUrl).cat.lineStream
          } catch {
            case exception: Exception =>
              throw new IllegalStateException(
                s"""|Verify that public repository https://hub.docker.com/r/$repositoryName exists.
                    |Either create the public repository including push access for the credentials in vault, or update
                    |`.withExecutableSettings()` in build.sbt adding either `buildDocker = false` or `pushDocker = false`.
                    |""".stripMargin,
                exception
              )
          }
        },
      )
    } else {
      List(
        DockerKeys.dockerPush := {
          ()
        },
      )
    }
  }

  // https://stackoverflow.com/questions/9819965/artifactory-snapshot-filename-handling
  private val buildTimestamp = System.currentTimeMillis() / 1000

  private def artifactoryResolver(isSnapshot: Boolean): Resolver = {
    val repoType = if (isSnapshot) "snapshot" else "release"
    val repoUrl =
      s"https://broadinstitute.jfrog.io/broadinstitute/libs-$repoType-local;build.timestamp=$buildTimestamp"
    val repoName = "artifactory-publish"
    repoName at repoUrl
  }

  private val artifactoryCredentialsFile =
    file("target/ci/resources/artifactory_credentials.properties").getAbsoluteFile

  private val artifactoryCredentials: Seq[Credentials] = {
    if (artifactoryCredentialsFile.exists)
      List(Credentials(artifactoryCredentialsFile))
    else
      Nil
  }

  val verifyArtifactoryCredentialsExist = taskKey[Unit]("Verify that the artifactory credentials file exists.")

  val artifactorySettings: Seq[Setting[_]] = List(
    publishTo := Option(artifactoryResolver(isSnapshot.value)),
    credentials ++= artifactoryCredentials,
  )

  val rootArtifactorySettings: Seq[Setting[_]] = List(
    verifyArtifactoryCredentialsExist := {
      if (!artifactoryCredentialsFile.exists) {
        throw new FileNotFoundException(
          s"""|${artifactoryCredentialsFile.toString}
              |The artifactory credentials file was not found.
              |Possibly the file was not rendered from vault,
              |or an `sbt clean` removed the /target directory.
              |""".stripMargin
        )
      }
    },
  )
}
