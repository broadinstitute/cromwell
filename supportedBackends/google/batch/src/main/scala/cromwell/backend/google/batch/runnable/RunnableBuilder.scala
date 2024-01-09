package cromwell.backend.google.batch.runnable

import com.google.cloud.batch.v1.Runnable.Container
import com.google.cloud.batch.v1.{Environment, Runnable, Volume}
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.batch.models.{BatchParameter, GcpBatchInput, GcpBatchOutput}
import cromwell.core.path.Path
import mouse.all.anySyntaxMouse

import scala.concurrent.duration.{Duration, DurationInt, FiniteDuration}
import scala.jdk.CollectionConverters._

/**
 * Utility singleton to create high level batch runnables.
 */
object RunnableBuilder {

  import RunnableLabels._
  import RunnableUtils._

  implicit class EnhancedRunnableBuilder(val builder: Runnable.Builder) extends AnyVal {

    /**
      * Only for use with docker images KNOWN to not have entrypoints already set,
      * or used with accompanying call to setEntrypoint("non-empty-string").
      *
      * Otherwise use the withEntrypointCommand() workaround below since the google issue listed in BA-6406 is not being
      * fixed.
      */
    def withCommand(command: String*): Runnable.Builder = {
      val container = builder.getContainerBuilder.addAllCommands(command.toList.asJava)
      builder.setContainer(container)
    }

    def withEntrypointCommand(command: String*): Runnable.Builder =
      builder
        .setContainer(
          builder.getContainerBuilder
            .setEntrypoint(
              command.headOption.getOrElse("")
            ) // set to blank string instead of null because batch does not support null
            .addAllCommands(
              command.drop(1).asJava
            )
        )

    def withFlags(flags: List[RunnableFlag]): Runnable.Builder =
      flags.foldLeft(builder) {
        case (acc, RunnableFlag.IgnoreExitStatus) => acc.setIgnoreExitStatus(true)
        case (acc, RunnableFlag.RunInBackground) => acc.setBackground(true)
        case (acc, RunnableFlag.AlwaysRun) => acc.setAlwaysRun(true)
      }

    def withEnvironment(environment: Map[String, String]): Runnable.Builder = {
      val env = Environment.newBuilder.putAllVariables(environment.asJava)
      builder.setEnvironment(env)
    }

    def withVolumes(volumes: List[Volume]): Runnable.Builder = {
      val formattedVolumes = volumes.map { volume =>
        val mountPath = volume.getMountPath

        val mountOptions = Option(volume.getMountOptionsList)
          .map(_.asScala)
          .filter(_.nonEmpty)
          .map(_.mkString(":", ",", ""))
          .getOrElse("")

        s"$mountPath:$mountPath$mountOptions"
      }

      builder.setContainer(
        builder.getContainerBuilder.addAllVolumes(formattedVolumes.asJava)
      )
    }

    def withLabels(labels: Map[String, String]): Runnable.Builder = builder.putAllLabels(labels.asJava)

    def withTimeout(timeout: Duration): Runnable.Builder = timeout match {
      case _: FiniteDuration =>
        builder.setTimeout(
          com.google.protobuf.Duration.newBuilder().setSeconds(timeout.toSeconds).build()
        )

      case _ => builder
    }

    def withAlwaysRun(alwaysRun: Boolean): Runnable.Builder = builder.setAlwaysRun(alwaysRun)

    def withRunInBackground(runInBackground: Boolean): Runnable.Builder = builder.setBackground(runInBackground)

    def scalaLabels: Map[String, String] = {
      val list = for {
        keyValueList <- Option(builder.getLabelsMap).toList
        keyValue <- keyValueList.asScala
      } yield keyValue
      list.toMap
    }
  }

  def withImage(image: String): Runnable.Builder =
    Runnable
      .newBuilder()
      .setContainer(Container.newBuilder.setImageUri(image))

  private def cloudSdkContainerBuilder: Container.Builder =
    Container.newBuilder.setImageUri(CloudSdkImage)

  def monitoringImageScriptRunnable(cloudPath: Path, containerPath: Path, volumes: List[Volume])(implicit
    gcsTransferConfiguration: GcsTransferConfiguration
  ): Runnable.Builder = {
    val command = RunnableCommands.localizeFile(cloudPath, containerPath)
    val labels = Map(Key.Tag -> Value.Localization)
    cloudSdkShellRunnable(command)(volumes = volumes, flags = List.empty, labels = labels)
  }

  def backgroundRunnable(image: String,
                         command: List[String],
                         environment: Map[String, String],
                         volumes: List[Volume]
  ): Runnable.Builder =
    withImage(image)
      .withEntrypointCommand(command: _*)
      .withRunInBackground(true)
      .withVolumes(volumes)
      .withEnvironment(environment)
      .withFlags(List(RunnableFlag.RunInBackground, RunnableFlag.IgnoreExitStatus))
      .withLabels(Map(Key.Tag -> Value.Monitoring))

  def terminateBackgroundRunnablesRunnable(): Runnable.Builder =
    cloudSdkShellRunnable(terminateAllBackgroundRunnablesCommand)(
      volumes = List.empty,
      flags = List(RunnableFlag.AlwaysRun),
      labels = Map(Key.Tag -> Value.Monitoring)
    )

  def gcsFileDeletionRunnable(cloudPath: String, volumes: List[Volume]): Runnable.Builder =
    cloudSdkShellRunnable(
      s"""gsutil rm '$cloudPath'"""
    )(
      volumes = volumes,
      flags = List(RunnableFlag.IgnoreExitStatus),
      labels = Map(Key.Tag -> Value.Monitoring)
    )

  def userRunnable(docker: String,
                   scriptContainerPath: String,
                   jobShell: String,
                   volumes: List[Volume],
                   dockerhubCredentials: (String, String)
  ): Runnable.Builder = {

    val container = (dockerhubCredentials._1, dockerhubCredentials._2) match {
      case (username, password) if username.nonEmpty && password.nonEmpty =>
        Container.newBuilder
          .setImageUri(docker)
          .setEntrypoint(jobShell)
          .addCommands(scriptContainerPath)
          .setUsername(username)
          .setPassword(password)
      case _ =>
        Container.newBuilder
          .setImageUri(docker)
          .setEntrypoint(jobShell)
          .addCommands(scriptContainerPath)
    }
    Runnable
      .newBuilder()
      .setContainer(container)
      .withVolumes(volumes)
      .putLabels(Key.Tag, Value.UserRunnable)
  }

  def checkForMemoryRetryRunnable(retryLookupKeys: List[String], volumes: List[Volume]): Runnable.Builder =
    cloudSdkShellRunnable(RunnableCommands.checkIfStderrContainsRetryKeys(retryLookupKeys))(
      volumes = volumes,
      flags = List(RunnableFlag.AlwaysRun),
      labels = Map(Key.Tag -> Value.RetryWithMoreMemory)
    ).withAlwaysRun(true)

  // Creates a Runnable that logs the docker command for the passed in runnable.
  def describeDocker(description: String, runnable: Runnable.Builder): Runnable.Builder =
    logTimestampedRunnable(
      s"Running $description: ${toDockerRun(runnable)}",
      List.empty,
      List.empty,
      runnable.scalaLabels
    )

  private def timestampedMessage(message: String): String =
    s"""printf '%s %s\\n' "$$(date -u '+%Y/%m/%d %H:%M:%S')" ${shellEscaped(message)}"""

  private def logTimestampedRunnable(message: String,
                                     volumes: List[Volume],
                                     flags: List[RunnableFlag],
                                     runnableLabels: Map[String, String]
  ): Runnable.Builder =
    // Uses the cloudSdk image as that image will be used for other operations as well.
    cloudSdkShellRunnable(
      timestampedMessage(message)
    )(volumes,
      flags,
      labels = runnableLabels collect {
        case (key, value) if key == Key.Tag => Key.Logging -> value
        case (key, value) => key -> value
      }
    ).withTimeout(timeout = 300.seconds)

  def cloudSdkRunnable: Runnable.Builder = Runnable.newBuilder.setContainer(cloudSdkContainerBuilder)

  def cloudSdkShellRunnable(shellCommand: String)(volumes: List[Volume],
                                                  flags: List[RunnableFlag],
                                                  labels: Map[String, String],
                                                  timeout: Duration = Duration.Inf
  ): Runnable.Builder =
    Runnable.newBuilder
      .setContainer(cloudSdkContainerBuilder)
      .withVolumes(volumes)
      .withLabels(labels)
      .withEntrypointCommand(
        "/bin/sh",
        "-c",
        if (shellCommand.contains("\n")) shellCommand |> RunnableCommands.multiLineCommand else shellCommand
      )
      .withFlags(flags)
      .withTimeout(timeout)

  def annotateTimestampedRunnable(description: String,
                                  loggingLabelValue: String,
                                  isAlwaysRun: Boolean = false
  )(volumes: List[Volume], runnables: List[Runnable.Builder]): List[Runnable.Builder] = {

    val flags = if (isAlwaysRun) List(RunnableFlag.AlwaysRun) else List()
    val labels = Map(Key.Logging -> loggingLabelValue)
    val starting = logTimestampedRunnable(s"Starting $description.", volumes, flags, labels)
    val done = logTimestampedRunnable(s"Done $description.", volumes, flags, labels)
    List(starting) ++ runnables ++ List(done)
  }

  /**
    * Returns a set of labels for a parameter.
    *
    * @param parameter Input or output parameter to label.
    * @return The labels.
    */
  def parameterLabels(parameter: BatchParameter): Map[String, String] =
    parameter match {
      case _: GcpBatchInput =>
        Map(
          Key.Tag -> Value.Localization,
          Key.InputName -> parameter.name
        )
      case _: GcpBatchOutput =>
        Map(
          Key.Tag -> Value.Delocalization,
          Key.OutputName -> parameter.name
        )
    }

  /** Creates a Runnable that describes the parameter localization or delocalization. */
  def describeParameter(parameter: BatchParameter,
                        volumes: List[Volume],
                        labels: Map[String, String]
  ): Runnable.Builder =
    parameter match {
      case _: GcpBatchInput =>
        val message = "Localizing input %s -> %s".format(
          shellEscaped(parameter.cloudPath),
          shellEscaped(parameter.containerPath)
        )
        logTimestampedRunnable(message, volumes, List.empty, labels)
      case _: GcpBatchOutput =>
        val message = "Delocalizing output %s -> %s".format(
          shellEscaped(parameter.containerPath),
          shellEscaped(parameter.cloudPath)
        )
        logTimestampedRunnable(message, volumes, List(RunnableFlag.AlwaysRun), labels)
    }

  // Converts an Runnable to a `docker run ...` command runnable in the shell.
  private[runnable] def toDockerRun(runnable: Runnable.Builder): String = {
    runnable.getContainer.getCommandsList.asScala.toList
      .map(cmd => shellEscaped(cmd))
      .mkString(" ")

    val commandArgs: String = Option(runnable.getContainerBuilder.getCommandsList) match {
      case Some(commands) =>
        commands.asScala map {
          case command if Option(command).isDefined => s" ${shellEscaped(command)}"
          case _ => ""
        } mkString ""
      case None => ""
    }

    val entrypointArg: String = Option(runnable.getContainerBuilder.getEntrypoint).filter(_.nonEmpty) match {
      case Some(entrypoint) => s" --entrypoint=${shellEscaped(entrypoint)}"
      case None => ""
    }

    val imageArg: String = Option(runnable.getContainerBuilder.getImageUri) match {
      case None => " <no docker image specified>"
      case Some(imageUri) => s" ${shellEscaped(imageUri)}"
    }

    val mountArgs: String = Option(runnable.getContainerBuilder.getVolumesList) match {
      case None => ""
      case Some(volumes) =>
        volumes.asScala map {
          case volume if Option(volume).isEmpty => ""
          case volume => s" -v ${shellEscaped(volume).replaceAll(":r[o|w]", "")}"
        } mkString ""
    }

    List("docker run", mountArgs, entrypointArg, imageArg, commandArgs).mkString
  }
}
