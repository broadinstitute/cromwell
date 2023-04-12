package cromwell.backend.google.pipelines.batch.runnable

import com.google.cloud.batch.v1.Runnable
import com.google.cloud.batch.v1.Runnable.Container
import cromwell.backend.google.pipelines.common.action.ActionUtils._

import scala.jdk.CollectionConverters._

/**
 * Utility singleton to create high level batch runnables.
 */
object RunnableBuilder {

  import RunnableLabels._

  implicit class EnhancedRunnableBuilder(val builder: Runnable.Builder) extends AnyVal {
    def withEntrypointCommand(command: String*): Runnable.Builder = {
      builder
        .setContainer(
          builder.getContainerBuilder
            .setEntrypoint(command.headOption.orNull)
            .addAllCommands(
              command.drop(1).asJava
            )
        )
    }

    def withAlwaysRun(alwaysRun: Boolean): Runnable.Builder = builder.setAlwaysRun(alwaysRun)

    def withRunInBackground(runInBackground: Boolean): Runnable.Builder = builder.setBackground(runInBackground)

    //  Runnable has labels in alpha.  Batch team adding to V1
//    def scalaLabels: Map[String, String] = {
//      val list = for {
//        keyValueList <- Option(runnable.getLabels).toList
//        keyValue <- keyValueList.asScala
//      } yield keyValue
//      list.toMap
//    }

//    def withVolumes(volumes: Volume): Runnable.Builder = runnable.toBuilder
//      .setContainer(Container.newBuilder.setVolumes(volumes))

    def withImage(image: String): Runnable.Builder = builder
      .setContainer(Container.newBuilder.setImageUri(image))
  }

  private def cloudSdkContainerBuilder: Container.Builder = Container.newBuilder.setImageUri(CloudSdkImage)

  //privateDockerKeyAndToken: Option[CreatePipelineDockerKeyAndToken],
  //fuseEnabled: Boolean)
  def userRunnable(docker: String,
                   command: String,
                   jobShell: String): Runnable.Builder = {

//    val dockerImageIdentifier = DockerImageIdentifier.fromString(docker)
//    val secret = for {
//      imageId <- dockerImageIdentifier.toOption
//      if DockerHub.isValidDockerHubHost(imageId.host) // This token only works for Docker Hub and not other repositories.
//      keyAndToken <- privateDockerKeyAndToken
//      s = new Secret().setKeyName(keyAndToken.key).setCipherText(keyAndToken.encryptedToken)
//    } yield s

    val container = Container.newBuilder
      .setImageUri(docker)
      .setEntrypoint(jobShell)
      .addCommands("-c") // TODO: Verify whether this is still required
      .addCommands(command)
    Runnable.newBuilder().setContainer(container)
    //.withLabels(labels)
    //.withTimeout(timeout)
  }


  //  Needs label support
  // Creates a Runnable that logs the docker command for the passed in action.
  def describeDocker(description: String, runnable: Runnable): Runnable = {
    logTimestampedAction(
      s"Running $description: ${RunnableBuilder.toDockerRun(runnable)}",
      Map.empty
    ).build()
  }

  private def timestampedMessage(message: String): String =
    s"""printf '%s %s\\n' "$$(date -u '+%Y/%m/%d %H:%M:%S')" ${shellEscaped(message)}"""

  private def logTimestampedAction(message: String,
                                   labels: Map[String, String]): Runnable.Builder = {
    // Uses the cloudSdk image as that image will be used for other operations as well.
    cloudSdkShellAction(
      timestampedMessage(message)
    )(labels)
  }

  // TODO: Use labels
  def cloudSdkShellAction(shellCommand: String)(labels: Map[String, String]): Runnable.Builder = {
    Runnable.newBuilder.setContainer(cloudSdkContainerBuilder)
      .withEntrypointCommand(
        "/bin/sh",
        "-c",
        shellCommand
      )
  }

    def annotateTimestampedRunnable(description: String, loggingLabelValue: String, isAlwaysRun: Boolean = false)
                                (actions: List[Runnable.Builder]): List[Runnable.Builder] = {

    val labels = Map(Key.Logging -> loggingLabelValue)
    val starting = logTimestampedAction(s"Starting $description.", labels).withAlwaysRun(isAlwaysRun)
    val done = logTimestampedAction(s"Done $description.", labels).withAlwaysRun(isAlwaysRun)
    List(starting) ++ actions ++ List(done)
  }

  // Converts an Runnable to a `docker run ...` command runnable in the shell.
  private[runnable] def toDockerRun(runnable: Runnable): String = {
    runnable.getContainer
      .getCommandsList
      .asScala
      .toList
      .map { cmd => shellEscaped(cmd) }
      .mkString(" ")
  }
}
