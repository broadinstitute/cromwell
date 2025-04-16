import Version.{cromwellVersion, Debug, Release, Snapshot, Standard}
import org.apache.ivy.Ivy
import org.apache.ivy.core.IvyPatternHelper
import org.apache.ivy.core.module.descriptor.{DefaultModuleDescriptor, MDArtifact}
import org.apache.ivy.plugins.resolver.IBiblioResolver
import sbt.Keys._
import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtdocker.DockerPlugin.autoImport._
import sbtdocker.Instruction

import java.io.FileNotFoundException
import scala.jdk.CollectionConverters._
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
      val tags = Version.buildType match {
        case Snapshot =>
          // Ordinary local build
          // Looks like `85-443a6fc-SNAP`
          Seq(version.value)
        case Release =>
          // Looks like `85`, `85-443a6fc`
          Seq(cromwellVersion, version.value)
        case Debug =>
          // Ordinary local build with debug stuff
          // Looks like `85-443a6fc-DEBUG`
          Seq(version.value)
        case Standard =>
          // Merge to `develop`
          // Looks like `85-443a6fc`, `latest`, `develop`
          // TODO: once we automate releases, `latest` should move to `Release`
          Seq(version.value, "latest", "develop")
      }

      val versionsCsv = tags.mkString(",")
      sys.env.getOrElse("CROMWELL_SBT_DOCKER_TAGS", versionsCsv).split(",")
    },
    docker / imageNames := dockerTags.value map { tag =>
      ImageName(namespace = Option("broadinstitute"), repository = name.value, tag = Option(tag))
    },
    docker / dockerfile := {
      // The assembly task generates a fat JAR file
      val artifact: File = assembly.value
      val artifactTargetPath = s"/app/${artifact.name}"
      val projectName = name.value
      val additionalDockerInstr: Seq[Instruction] = (dockerCustomSettings ?? Nil).value

      new Dockerfile {
        from("us.gcr.io/broad-dsp-gcr-public/base/jre:17-debian")
        expose(8000)
        add(artifact, artifactTargetPath)
        runRaw(s"ln -s $artifactTargetPath /app/$projectName.jar")

        // Extra tools in debug mode only
        if (Version.buildType == Debug) {
          addInstruction(installDebugFacilities(version.value))
        }

        // Add a custom java opt, this avoids the following error (from Akka):
        //     class com.typesafe.sslconfig.ssl.DefaultHostnameVerifier (in unnamed module @0x5594a1b5)
        //     cannot access class sun.security.util.HostnameChecker (in module java.base)
        //     because module java.base does not export sun.security.util to unnamed module @0x5594a1b5
        // See https://docs.oracle.com/en/java/javase/17/migrate/migrating-jdk-8-later-jdk-releases.html#GUID-2F61F3A9-0979-46A4-8B49-325BA0EE8B66
        // TODO remove this once we upgrade Akka past 2.5
        val addOpensJavaOpt = "--add-opens=java.base/sun.security.util=ALL-UNNAMED"

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
          s"java $${JAVA_OPTS} ${addOpensJavaOpt} -jar /app/$projectName.jar $${${projectName.toUpperCase.replaceAll("-", "_")}_ARGS} $${*}",
          "--"
        )
        // for each custom setting (instruction) run addInstruction()
        additionalDockerInstr.foreach(addInstruction)
      }
    },
    docker / buildOptions := BuildOptions(
      cache = false,
      removeIntermediateContainers = BuildOptions.Remove.Always
    )
  )

  /**
    * Install packages needed for debugging when shelled in to a running image.
    *
    * This includes:
    * - The JDK, which includes tools like `jstack` not present in the JRE
    * - The YourKit Java Profiler
    * - Various Linux system & development utilities
    *
    * @return Instruction to run in the build
    */
  def installDebugFacilities(displayVersion: String): Instruction = {
    import sbtdocker.Instructions
    import java.time.{ZoneId, ZonedDateTime}
    import java.time.format.DateTimeFormatter

    val buildTime =
      ZonedDateTime.now(ZoneId.systemDefault()).format(DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss z"))

    // It is optimal to use a single `Run` instruction to minimize the number of layers in the image.
    //
    // Documentation:
    // - https://www.yourkit.com/docs/java-profiler/2024.3/help/docker_broker.jsp#setup
    Instructions.Run(
      s"""apt-get update -qq && \\
         |apt-get install -qq --no-install-recommends file gpg htop jq less nload unzip vim && \\
         |rm -rf /var/lib/apt/lists/* && \\
         |mkdir /tmp/docker-build-cache && \\
         |wget -q https://www.yourkit.com/download/docker/YourKit-JavaProfiler-2024.3-docker.zip -P /tmp/docker-build-cache/ && \\
         |unzip /tmp/docker-build-cache/YourKit-JavaProfiler-2024.3-docker.zip -d /tmp/docker-build-cache && \\
         |mkdir -p /usr/local/YourKit-JavaProfiler-2024.3/bin/ && \\
         |cp -R /tmp/docker-build-cache/YourKit-JavaProfiler-2024.3/bin/linux-x86-64/ /usr/local/YourKit-JavaProfiler-2024.3/bin/linux-x86-64/ && \\
         |rm -rf /tmp/docker-build-cache && \\
         |echo "Version $displayVersion built at $buildTime" > /etc/motd && \\
         |echo "[ ! -z "\\$$TERM" -a -r /etc/motd ] && cat /etc/motd" > /etc/bash.bashrc
         |""".stripMargin
    )
  }

  def dockerPushSettings(pushEnabled: Boolean): Seq[Setting[_]] =
    if (pushEnabled) {
      List(
        dockerPushCheck := {
          val projectName = name.value
          val repositoryName = s"broadinstitute/$projectName"
          val repositoryUrl = s"https://registry.hub.docker.com/v2/repositories/$repositoryName/"
          try
            url(repositoryUrl).cat.lineStream
          catch {
            case exception: Exception =>
              throw new IllegalStateException(
                s"""|Verify that public repository https://hub.docker.com/r/$repositoryName exists.
                    |Either create the public repository including push access for the credentials in vault, or update
                    |`.withExecutableSettings()` in build.sbt adding either `buildDocker = false` or `pushDocker = false`.
                    |""".stripMargin,
                exception
              )
          }
        }
      )
    } else {
      List(
        DockerKeys.dockerPush := {
          Map.empty[sbtdocker.ImageName, sbtdocker.ImageDigest]
        }
      )
    }

  private val broadArtifactoryResolver: Resolver =
    "Broad Artifactory" at
      "https://broadinstitute.jfrog.io/broadinstitute/libs-release/"

  private val broadArtifactoryResolverSnap: Resolver =
    "Broad Artifactory Snapshots" at
      "https://broadinstitute.jfrog.io/broadinstitute/libs-snapshot-local/"

  // https://stackoverflow.com/questions/9819965/artifactory-snapshot-filename-handling
  private val buildTimestamp = System.currentTimeMillis() / 1000

  private val broadArtifactoryLocalResolver: Resolver =
    "Broad Artifactory Local" at
      s"https://broadinstitute.jfrog.io/broadinstitute/libs-release-local;build.timestamp=$buildTimestamp/"

  val additionalResolvers = List(
    broadArtifactoryResolver,
    broadArtifactoryResolverSnap
  ) ++ Resolver.sonatypeOssRepos("releases")

  private val artifactoryCredentialsFile =
    file("target/ci/resources/artifactory_credentials.properties").getAbsoluteFile

  private val artifactoryCredentials: Seq[Credentials] =
    if (artifactoryCredentialsFile.exists)
      List(Credentials(artifactoryCredentialsFile))
    else
      Nil

  // BT-250 Check if publishing will fail due to already published artifacts
  val checkAlreadyPublished = taskKey[Boolean]("Verifies if publishing has already occurred")
  val errorIfAlreadyPublished = taskKey[Unit]("Fails the build if publishing has already occurred")

  private case class CromwellMDArtifactType(artifactType: String,
                                            artifactExtension: String,
                                            classifierOption: Option[String]
  )

  /**
    * The types of MDArtifacts published by this sbt build.
    *
    * This static list is an alternative to retrieving the types from sbt's `publishConfiguration` and its `artifacts`.
    * That `publishConfiguration` unfortunately depends on `compile` so it takes several minutes to run.
    */
  private val cromwellMDArtifactTypes = List(
    CromwellMDArtifactType("pom", "pom", None),
    CromwellMDArtifactType("jar", "jar", None),
    CromwellMDArtifactType("src", "jar", Option("sources")),
    CromwellMDArtifactType("doc", "jar", Option("javadoc"))
  )

  /**
    * Retrieve the IBiblioResolver from sbt's Ivy setup.
    */
  private def getIBiblioResolver(ivy: Ivy): IBiblioResolver =
    ivy.getSettings.getResolver(broadArtifactoryResolver.name) match {
      case iBiblioResolver: IBiblioResolver => iBiblioResolver
      case other => sys.error(s"Expected an IBiblioResolver, got $other")
    }

  /**
    * Maps an sbt artifact to the Apache Ivy artifact type.
    */
  private def makeMDArtifact(moduleDescriptor: DefaultModuleDescriptor)(
    cromwellMDArtifactType: CromwellMDArtifactType
  ): MDArtifact =
    new MDArtifact(
      moduleDescriptor,
      moduleDescriptor.getModuleRevisionId.getName,
      cromwellMDArtifactType.artifactType,
      cromwellMDArtifactType.artifactExtension,
      null,
      cromwellMDArtifactType.classifierOption.map("classifier" -> _).toMap.asJava
    )

  /**
    * Returns true and prints out an error if an artifact already exists.
    */
  private def existsMDArtifact(resolver: IBiblioResolver, log: Logger)(mdArtifact: MDArtifact): Boolean = {
    val exists = resolver.exists(mdArtifact)
    if (exists) {
      val pattern = resolver.getRoot + resolver.getPattern
      val urlString = IvyPatternHelper.substitute(pattern, mdArtifact)
      log.error(s"Already published: $urlString")
    }
    exists
  }

  val publishingSettings: Seq[Setting[_]] = List(
    publishTo := Option(broadArtifactoryLocalResolver),
    credentials ++= artifactoryCredentials,
    checkAlreadyPublished := {
      val module = ivyModule.value
      val log = streams.value.log

      module.withModule(log) { case (ivy, moduleDescriptor, _) =>
        val resolver = getIBiblioResolver(ivy)
        cromwellMDArtifactTypes
          .map(makeMDArtifact(moduleDescriptor))
          .map(existsMDArtifact(resolver, log))
          .exists(identity)
      }
    },
    errorIfAlreadyPublished := {
      if (checkAlreadyPublished.value) {
        sys.error(
          s"Some ${version.value} artifacts were already published and will need to be manually deleted. " +
            "See the errors above for the list of published artifacts."
        )
      }
    }
  )

  val verifyArtifactoryCredentialsExist = taskKey[Unit]("Verify that the artifactory credentials file exists.")

  val rootPublishingSettings: Seq[Setting[_]] = List(
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
    }
  )
}
