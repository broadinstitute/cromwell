import java.io.FileNotFoundException

import Version.cromwellVersion
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtdocker.DockerPlugin.autoImport._
import ContinuousIntegration._
import sbtdocker.Instruction

import scala.sys.process._

object Publishing {

  val dockerTags = settingKey[Seq[String]]("The tags for docker builds.")
  val dockerPushCheck = taskKey[Unit]("Check that the public repository exists in DockerHub.")

  val dockerCustomSettings = settingKey[Seq[Instruction]]("Additional instructions for docker image.")

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
      val additionalDockerInstr: Seq[Instruction] = dockerCustomSettings.value

      new Dockerfile {
        from("us.gcr.io/broad-dsp-gcr-public/base/jre:8-debian")
        expose(8000)
        add(artifact, artifactTargetPath)
        runRaw(s"ln -s $artifactTargetPath /app/$projectName.jar")

        /*
        If you use the 'exec' form for an entry point, shell processing is not performed and
        environment variable substitution does not occur.  Thus we have to /bin/bash here
        and pass along any subsequent command line arguments
        See https://docs.docker.com/engine/reference/builder/#/entrypoint

        Notes and warnings on JAVA_OPTS in docker-compose YAML files:

        Setting JAVA_OPTS in a docker-compose YAML requires a combination of passing and parsing values through:
        1. docker-compose
        2. environment key/values
        3. bash parameter expansion without full bash command parsing

        Examples:
        - docker-compose splits the env key/value on the first `=`:
          - If in the docker-compose YAML:              my_key=my_val1=my_al2
          - Then in the environment key "my_key" value: my_val1=my_val2
        - For legibility the YAMLs use chomped-block-scalar for JAVA_OPTS: https://yaml-multiline.info/
          - Newlines are removed, and the yaml value is compressed into a space separated string.
          - Again, only the first `=` sign is used by docker-compose to separate the env key/value pair.
          - The newlines-that-are-now-spaces get included into the env values.
        - Do not quote args! The quotes will be passed into the env-value by docker-compose. Bash will not remove them.
          - If in the docker-compose YAML: -Dmy_key="my_val"
          - Then in the running program:   conf.getString("my_key") == "\"my_val\"")
        - Spaces are ok in YAML values, but NOT env values! The JAVA_OPTS in the `entrypoint` below does not "quote".
          - If the env key "ENV_PATH" has value: /path/dir with spaces/file.conf
          - If in the docker-compose YAML:       JAVA_OPTS=-Dconfig.file=${ENV_PATH} -Dfoo=bar
          - Then `entryPoint` bash runs command: "java" "-Dconf.file=/path/dir" "with" "spaces/file.conf" "-Dfoo=bar"
        - If needed use $$ to escape dollar signs: https://docs.docker.com/compose/compose-file/#variable-substitution
          - Don't have a current example where this is required
          - Often one can just allow docker-compose to pass or perform environment variable substitution
         */
        entryPoint(
          "/bin/bash",
          "-c",
          s"java $${JAVA_OPTS} -jar /app/$projectName.jar $${${projectName.toUpperCase.replaceAll("-", "_")}_ARGS} $${*}",
          "--"
        )
        // for each custom setting (instruction) run addInstruction()
        additionalDockerInstr.foreach(addInstruction)
      }
    },
    buildOptions in docker := BuildOptions(
      cache = false,
      removeIntermediateContainers = BuildOptions.Remove.Always
    ),
    dockerCustomSettings in ThisBuild := Nil, // setting the default value
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
          Map.empty[sbtdocker.ImageName,sbtdocker.ImageDigest]
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
