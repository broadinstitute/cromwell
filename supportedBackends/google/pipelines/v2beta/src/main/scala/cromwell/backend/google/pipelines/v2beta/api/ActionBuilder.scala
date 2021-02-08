package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount, Secret}
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.action.ActionLabels._
import cromwell.backend.google.pipelines.common.action.ActionUtils._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineDockerKeyAndToken
import cromwell.backend.google.pipelines.common.{PipelinesApiInput, PipelinesApiOutput, PipelinesParameter}
import cromwell.backend.google.pipelines.v2beta.LifeSciencesFactory
import cromwell.core.path.Path
import cromwell.docker.DockerImageIdentifier
import cromwell.docker.registryv2.flows.dockerhub.DockerHub
import mouse.all._

import scala.collection.JavaConverters._
import scala.concurrent.duration._

/**
  * Utility singleton to create high level actions.
  */
object ActionBuilder {
  implicit class EnhancedAction(val action: Action) extends AnyVal {
    /**
      * Only for use with docker images KNOWN to not have entrypoints already set,
      * or used with accompanying call to setEntrypoint("non-empty-string").
      *
      * Otherwise use the withEntrypointCommand() workaround below since the google issue listed in BA-6406 is not being
      * fixed.
      */
    def withCommand(command: String*): Action = action.setCommands(command.toList.asJava)

    /**
      * Useful for any externally provided images that _might_ have entrypoints already set. This is a workaround for
      * the issue detailed in BA-6406. See underlying google issue in that ticket for more info.
      */
    def withEntrypointCommand(command: String*): Action = {
      action
        .setEntrypoint(command.headOption.orNull)
        .setCommands(
          Option(command.drop(1))
            .filter(_.nonEmpty)
            .map(_.asJava)
            .orNull
        )
    }

    def withMounts(mounts: List[Mount]): Action = action.setMounts(mounts.asJava)
    def withLabels(labels: Map[String, String]): Action = action.setLabels(labels.asJava)
    def withTimeout(timeout: Duration): Action = timeout match {
      case fd: FiniteDuration => action.setTimeout(fd.toSeconds + "s")
      case _ => action
    }
    def withIgnoreExitStatus(ignoreExitStatus: Boolean): Action = action.setIgnoreExitStatus(ignoreExitStatus)
    def withRunInBackground(runInBackground: Boolean): Action = action.setRunInBackground(runInBackground)
    def withAlwaysRun(alwaysRun: Boolean): Action = action.setAlwaysRun(alwaysRun)
    def withEnableFuse(enableFuse: Boolean): Action = action.setEnableFuse(enableFuse)
    def withPublishExposedPorts(publishExposedPorts: Boolean): Action = action.setPublishExposedPorts(publishExposedPorts)
    def withDisableImagePrefetch(disableImagePrefetch: Boolean): Action = action.setDisableImagePrefetch(disableImagePrefetch)

    def scalaLabels: Map[String, String] = {
      val list = for {
        keyValueList <- Option(action.getLabels).toList
        keyValue <- keyValueList.asScala
      } yield keyValue
      list.toMap
    }
  }

  def cloudSdkAction: Action = new Action().setImageUri(LifeSciencesFactory.CloudSdkImage)

  def withImage(image: String): Action = new Action()
    .setImageUri(image)

  def monitoringImageScriptAction(cloudPath: Path, containerPath: Path, mounts: List[Mount])
                                 (implicit gcsTransferConfiguration: GcsTransferConfiguration): Action = {
    val command = ActionCommands.localizeFile(cloudPath, containerPath)
    val labels = Map(Key.Tag -> Value.Localization)
    ActionBuilder.cloudSdkShellAction(command)(mounts = mounts, labels = labels)
  }

  def backgroundAction(image: String,
                       command: List[String],
                       environment: Map[String, String],
                       mounts: List[Mount],
                      ): Action = {
    new Action()
      .setImageUri(image)
      .withEntrypointCommand(command: _*)
      .withRunInBackground(true)
      .withIgnoreExitStatus(true)
      .withMounts(mounts)
      .setEnvironment(environment.asJava)
      .withLabels(Map(Key.Tag -> Value.Monitoring))
      .setPidNamespace(backgroundActionPidNamespace)
  }

  def terminateBackgroundActionsAction(): Action =
    cloudSdkShellAction(terminateAllBackgroundActionsCommand)(labels = Map(Key.Tag -> Value.Monitoring))
      .withAlwaysRun(true)
      .setPidNamespace(backgroundActionPidNamespace)

  def gcsFileDeletionAction(cloudPath: String): Action =
    cloudSdkShellAction(
      s"""gsutil rm '$cloudPath'"""
    )(labels = Map(Key.Tag -> Value.Monitoring))
      .withIgnoreExitStatus(true)

  def userAction(docker: String,
                 scriptContainerPath: String,
                 mounts: List[Mount],
                 jobShell: String,
                 privateDockerKeyAndToken: Option[CreatePipelineDockerKeyAndToken],
                 fuseEnabled: Boolean): Action = {

    val dockerImageIdentifier = DockerImageIdentifier.fromString(docker)

    val secret = for {
      imageId <- dockerImageIdentifier.toOption
      if DockerHub.isValidDockerHubHost(imageId.host) // This token only works for Docker Hub and not other repositories.
      keyAndToken <- privateDockerKeyAndToken
      s = new Secret().setKeyName(keyAndToken.key).setCipherText(keyAndToken.encryptedToken)
    } yield s

    new Action()
      .setImageUri(docker)
      .setCommands(List(scriptContainerPath).asJava)
      .setMounts(mounts.asJava)
      .setEntrypoint(jobShell)
      .setLabels(Map(Key.Tag -> Value.UserAction).asJava)
      .setCredentials(secret.orNull)
      .setEnableFuse(fuseEnabled)
  }

  def checkForMemoryRetryAction(retryLookupKeys: List[String], mounts: List[Mount]): Action = {
    cloudSdkShellAction(ActionCommands.checkIfStderrContainsRetryKeys(retryLookupKeys))(
      mounts = mounts,
      labels = Map(Key.Tag -> Value.RetryWithMoreMemory)
    ).withAlwaysRun(true)
  }

  def cloudSdkShellAction(shellCommand: String)(mounts: List[Mount] = List.empty,
                                                labels: Map[String, String] = Map.empty,
                                                timeout: Duration = Duration.Inf): Action =
    cloudSdkAction
      .withEntrypointCommand(
        "/bin/sh",
        "-c",
        if (shellCommand.contains("\n")) shellCommand |> ActionCommands.multiLineCommand else shellCommand
      )
      .withMounts(mounts)
      .withLabels(labels)
      .withTimeout(timeout)

  /**
    * Returns a set of labels for a parameter.
    *
    * @param pipelinesParameter Input or output parameter to label.
    * @return The labels.
    */
  def parameterLabels(pipelinesParameter: PipelinesParameter): Map[String, String] = {
    pipelinesParameter match {
      case _: PipelinesApiInput =>
        Map(
          Key.Tag -> Value.Localization,
          Key.InputName -> pipelinesParameter.name
        )
      case _: PipelinesApiOutput =>
        Map(
          Key.Tag -> Value.Delocalization,
          Key.OutputName -> pipelinesParameter.name
        )
    }
  }

  /**
    * Surrounds the list of Actions with a pair of starting and done Actions.
    *
    * @param description       Description of the list of Actions.
    * @param loggingLabelValue An entry from Value that describes the list of Actions.
    * @param isAlwaysRun       If true the pair of starting and done Actions will always be run.
    * @param actions           The list of Actions to surround.
    * @return The starting Action, the passed in list, and then a done Action.
    */
  def annotateTimestampedActions(description: String, loggingLabelValue: String, isAlwaysRun: Boolean = false)
                                (actions: List[Action]): List[Action] = {
    val labels = Map(Key.Logging -> loggingLabelValue)
    val starting = List(logTimestampedAction(s"Starting $description.", labels).withAlwaysRun(isAlwaysRun))
    val done = List(logTimestampedAction(s"Done $description.", labels).withAlwaysRun(isAlwaysRun))
    starting ++ actions ++ done
  }

  /** Creates an Action that describes the parameter localization or delocalization. */
  def describeParameter(pipelinesParameter: PipelinesParameter,
                        actionLabels: Map[String, String]): Action = {
    pipelinesParameter match {
      case _: PipelinesApiInput =>
        val message = "Localizing input %s -> %s".format(
          shellEscaped(pipelinesParameter.cloudPath),
          shellEscaped(pipelinesParameter.containerPath),
        )
        ActionBuilder.logTimestampedAction(message, actionLabels)
      case _: PipelinesApiOutput =>
        val message = "Delocalizing output %s -> %s".format(
          shellEscaped(pipelinesParameter.containerPath),
          shellEscaped(pipelinesParameter.cloudPath),
        )
        ActionBuilder.logTimestampedAction(message, actionLabels).withAlwaysRun(true)
    }
  }

  /** Creates an Action that logs the docker command for the passed in action. */
  def describeDocker(description: String, action: Action): Action = {
    ActionBuilder.logTimestampedAction(
      s"Running $description: ${ActionBuilder.toDockerRun(action)}",
      action.scalaLabels
    )
  }

  def timestampedMessage(message: String): String =
    s"""printf '%s %s\\n' "$$(date -u '+%Y/%m/%d %H:%M:%S')" ${shellEscaped(message)}"""

  /**
    * Creates an Action that logs the time as UTC plus prints the message. The original actionLabels will also be
    * applied to the logged action, except that Key.Tag -> some-value will be replaced with Key.Logging -> some-value.
    *
    * Note that these log actions have a timeout of 300 seconds. That's obviously a huge timeout, and is more of a
    * safety net in case these actions get stuck.
    *
    * @param message      Message to output.
    * @param actionLabels Labels from the original Action to modify and apply to the logging Action.
    * @return A new Action that will log the time and print the message.
    */
  private def logTimestampedAction(message: String,
                                   actionLabels: Map[String, String]): Action = {
    // Uses the cloudSdk image as that image will be used for other operations as well.
    cloudSdkShellAction(
      timestampedMessage(message)
    )(
      labels = actionLabels collect {
        case (key, value) if key == Key.Tag => Key.Logging -> value
        case (key, value) => key -> value
      },
      timeout = 300.seconds
    )
  }

  /** Converts an Action to a `docker run ...` command runnable in the shell. */
  private[api] def toDockerRun(action: Action): String = {
    val commandArgs: String = Option(action.getCommands) match {
      case Some(commands) =>
        commands.asScala map {
          case command if Option(command).isDefined => s" ${shellEscaped(command)}"
          case _ => ""
        } mkString ""
      case None => ""
    }

    val entrypointArg: String = Option(action.getEntrypoint) match {
      case Some(entrypoint) => s" --entrypoint=${shellEscaped(entrypoint)}"
      case None => ""
    }

    val environmentArgs: String = Option(action.getEnvironment) match {
      case Some(environment) =>
        environment.asScala map {
          case (key, value) if Option(key).isDefined && Option(value).isDefined => s" -e ${shellEscaped(s"$key:$value")}"
          case (key, _) if Option(key).isDefined => s" -e ${shellEscaped(key)}"
          case _ => ""
        } mkString ""
      case None => ""
    }

    val imageArg: String = Option(action.getImageUri) match {
      case None => " <no docker image specified>"
      case Some(imageUri) => s" ${shellEscaped(imageUri)}"
    }

    val mountArgs: String = Option(action.getMounts) match {
      case None => ""
      case Some(mounts) =>
        mounts.asScala map {
          case mount if Option(mount).isEmpty => ""
          case mount if mount.getReadOnly => s" -v ${shellEscaped(s"/mnt/${mount.getDisk}:${mount.getPath}:ro")}"
          case mount => s" -v ${shellEscaped(s"/mnt/${mount.getDisk}:${mount.getPath}")}"
        } mkString ""
    }

    val nameArg: String = Option(action.getContainerName) match {
      case None => ""
      case Some(name) => s" --name ${shellEscaped(name)}"
    }

    val pidNamespaceArg: String = Option(action.getPidNamespace) match {
      case Some(pidNamespace) => s" --pid=${shellEscaped(pidNamespace)}"
      case None => ""
    }

    val portMappingArgs: String = Option(action.getPortMappings) match {
      case Some(portMappings) =>
        portMappings.asScala map {
          case (key, value) if Option(key).isDefined => s" -p ${shellEscaped(s"$key:$value")}"
          case (_, value) => s" -p $value"
        } mkString ""
      case None => ""
    }

    val publishExposedPortsFlag: String = if (action.getPublishExposedPorts) " -P" else ""

    Array("docker run",
      nameArg,
      mountArgs,
      environmentArgs,
      pidNamespaceArg,
      publishExposedPortsFlag,
      portMappingArgs,
      entrypointArg,
      imageArg,
      commandArgs,
    ).mkString
  }
}
