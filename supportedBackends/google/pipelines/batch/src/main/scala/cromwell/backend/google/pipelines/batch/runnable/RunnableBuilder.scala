package cromwell.backend.google.pipelines.batch.runnable


//import com.google.cloud.batch.v1.{Runnable, Volume}
//import com.google.cloud.batch.v1.Runnable.Container
//import scala.concurrent.duration._
//import cromwell.backend.google.pipelines.common.action.ActionCommands
//import cromwell.backend.google.pipelines.common.action.ActionUtils._
//import cromwell.backend.google.pipelines.common.action.ActionLabels._

//import mouse.all._

//import scala.jdk.CollectionConverters._
//import scala.concurrent.duration._

/**
 * Utility singleton to create high level batch runnables.
 */

object RunnableBuilder {
/*
  implicit class EnhancedRunnable(val runnable: Runnable) extends AnyVal {


    def withCommand(commands: List[String]): Runnable.Builder = {


      new Runnable()
      .toBuilder
        .setContainer(Container.newBuilder.addAllCommands(commands.asJava))
    }
    /*
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
    }*/

    //def withRunInBackground(runInBackground: Boolean): Runnable = runnable.setBackground(runInBackground)

    def withAlwaysRun(alwaysRun: Boolean): Runnable = runnable.withAlwaysRun(alwaysRun=true)


    /*  Runnable has labels in alpha.  Batch team adding to V1
    def scalaLabels: Map[String, String] = {
      val list = for {
        keyValueList <- Option(runnable.getLabels).toList
        keyValue <- keyValueList.asScala
      } yield keyValue
      list.toMap
    }

     */

    def withVolumes(volumes: Volume): Runnable.Builder = new Runnable()
      .toBuilder
      .setContainer(Container.newBuilder.setVolumes(volumes))

  }

  def cloudSdkAction: Runnable.Builder = new Runnable()
    .toBuilder
    .setContainer(Container.newBuilder.setImageUri(CloudSdkImage))


  def withImage(image: String): Runnable.Builder = new Runnable()
    .toBuilder
    .setContainer(Container.newBuilder.setImageUri(image))




  def userRunnable(docker: String,
                 scriptContainerPath: String,
                   volumes: List[Volume],
                 jobShell: String)
                 //privateDockerKeyAndToken: Option[CreatePipelineDockerKeyAndToken],
                 //fuseEnabled: Boolean)
  : Runnable.Builder = {

    //val dockerImageIdentifier = DockerImageIdentifier.fromString(docker)

    /*
    val secret = for {
      imageId <- dockerImageIdentifier.toOption
      if DockerHub.isValidDockerHubHost(imageId.host) // This token only works for Docker Hub and not other repositories.
      keyAndToken <- privateDockerKeyAndToken
      s = new Secret().setKeyName(keyAndToken.key).setCipherText(keyAndToken.encryptedToken)
    } yield s
    */


    new Runnable().toBuilder.setContainer(Container.newBuilder.setImageUri(docker).addCommands(scriptContainerPath).setEntrypoint(jobShell))
  }

  /*
  def cloudSdkShellAction(shellCommand: String)(volumes: Volume,
                                                labels: Map[String, String] = Map.empty,
                                                timeout: Duration = Duration.Inf): Runnable.Builder = {
    //withVolumes(volumes)
    cloudSdkAction
      .setContainer(Container.newBuilder.addVolumes(volumes))
      //.withEntrypointCommand(
      //  "/bin/sh",
      //  "-c",
      //  if (shellCommand.contains("\n")) shellCommand |> ActionCommands.multiLineCommand else shellCommand
      //)
      //.withVolumes(volumes)
  }
  //.withLabels(labels)
      //.withTimeout(timeout)

*/

  /*  Needs label support
  /** Creates an Action that logs the docker command for the passed in action. */
  def describeDocker(description: String, runnable: Runnable): Runnable = {
    RunnableBuilder.logTimestampedAction(
      s"Running $description: ${RunnableBuilder.toDockerRun(runnable)}",
      runnable.scalaLabels
    )
  }

   */

  def timestampedMessage(message: String): String =
    s"""printf '%s %s\\n' "$$(date -u '+%Y/%m/%d %H:%M:%S')" ${shellEscaped(message)}"""

  private def logTimestampedAction(message: String,
                                   actionLabels: Map[String, String]): Runnable.Builder = {
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
    val starting = List(logTimestampedAction(s"Starting $description.", labels).buildPartial.withAlwaysRun(isAlwaysRun))
    val done = List(logTimestampedAction(s"Done $description.", labels).buildPartial.withAlwaysRun(isAlwaysRun))
    starting ++ actions ++ done
  }


  /*
  /** Converts an Action to a `docker run ...` command runnable in the shell. */
  private[api] def toDockerRun(runnable: Runnable): String = {
    val commandArgs: String = Option(runnable.getCommands) match {
      case Some(commands) =>
        commands.asScala map {
          case command if Option(command).isDefined => s" ${shellEscaped(command)}"
          case _ => ""
        } mkString ""
      case None => ""
    }

   */

 */
}
