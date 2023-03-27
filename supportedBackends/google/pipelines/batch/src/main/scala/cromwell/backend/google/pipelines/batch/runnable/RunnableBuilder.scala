package cromwell.backend.google.pipelines.batch.runnable


import com.google.cloud.batch.v1.Runnable
import com.google.cloud.batch.v1.Runnable.Container
import scala.concurrent.duration._
import cromwell.backend.google.pipelines.common.action.ActionCommands
import cromwell.backend.google.pipelines.common.action.ActionUtils._
import cromwell.backend.google.pipelines.common.action.ActionLabels._
import mouse.all._
import scala.jdk.CollectionConverters._

/**
 * Utility singleton to create high level batch runnables.
 */

object RunnableBuilder {

  implicit class EnhancedRunnable(val runnable: Runnable) extends AnyVal {


    def withEntrypointCommand(command: String*): Runnable = {
      Runnable
        .newBuilder
        .setContainer(Container
          .newBuilder
          .setEntrypoint(command.headOption.orNull)
          .setCommands(
            Option(command.drop(1))
              .filter(_.nonEmpty)
              .map(_.asJava)
              .orNull
          )
        )
    }


    //def withRunInBackground(runInBackground: Boolean): Runnable = runnable.setBackground(runInBackground)
    def withAlwaysRun(alwaysRun: Boolean): Runnable = runnable.withAlwaysRun(alwaysRun=true)

  }

  def cloudSdkAction: Runnable = Runnable
    .newBuilder
    .setContainer((Container
      .newBuilder
      .setImageUri(CloudSdkImage)
      .setEntrypoint("/bin/sh")
      .addCommands("-c")
      .addCommands("echo Hello World!")
      .build
      ))
    .build

  def cloudSdkShellAction(shellCommand: String)(mounts: List[Mount] = List.empty,
                                                labels: Map[String, String] = Map.empty,
                                                timeout: Duration = Duration.Inf): Runnable =
    cloudSdkAction
      .withEntrypointCommand(
        "/bin/sh",
        "-c",
        if (shellCommand.contains("\n")) shellCommand |> ActionCommands.multiLineCommand else shellCommand
      )
      .withMounts(mounts)
      .withLabels(labels)
      .withTimeout(timeout)

  def timestampedMessage(message: String): String =
    s"""printf '%s %s\\n' "$$(date -u '+%Y/%m/%d %H:%M:%S')" ${shellEscaped(message)}"""

  private def logTimestampedAction(message: String,
                                   actionLabels: Map[String, String]): Runnable = {
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

  def annotateTimestampedActions(description: String, loggingLabelValue: String, isAlwaysRun: Boolean = false)
                                (actions: List[Runnable]): List[Runnable] = {
    val labels = Map(Key.Logging -> loggingLabelValue)
    val starting = List(logTimestampedAction(s"Starting $description.", labels).withAlwaysRun(isAlwaysRun))
    val done = List(logTimestampedAction(s"Done $description.", labels).withAlwaysRun(isAlwaysRun))
    starting ++ actions ++ done
  }

}
