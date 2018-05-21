package cromwell.backend.google.pipelines.v2alpha1.api

import akka.http.scaladsl.model.ContentTypes
import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.PipelinesApiFileOutput
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Gsutil.ContentTypeTextHeader
import cromwell.backend.google.pipelines.v2alpha1.api.ActionFlag.ActionFlag
import mouse.all._
import org.apache.commons.text.StringEscapeUtils.ESCAPE_XSI

import scala.collection.JavaConverters._

/**
  * Utility singleton to create high level actions.
  */
object ActionBuilder {
  implicit class EnhancedAction(val action: Action) extends AnyVal {
    private def javaFlags(flags: List[ActionFlag]) = flags.map(_.toString).asJava

    def withCommand(command: String*): Action = action.setCommands(command.toList.asJava)
    def withFlags(flags: List[ActionFlag]): Action = action.setFlags(flags |> javaFlags)
    def withMounts(mounts: List[Mount]): Action = action.setMounts(mounts.asJava)
    def withLabels(labels: Map[String, String]): Action = action.setLabels(labels.asJava)
  }
  
  object Gsutil {
    private val contentTypeText = ContentTypes.`text/plain(UTF-8)`.toString()
    val ContentTypeTextHeader = s"Content-Type: $contentTypeText"
  }

  private val cloudSdkImage = "google/cloud-sdk:alpine"
  def cloudSdkAction: Action = new Action().setImageUri(cloudSdkImage)
  
  def withImage(image: String) = new Action()
    .setImageUri(image)

  def userAction(docker: String, scriptContainerPath: String, mounts: List[Mount]): Action = {
    new Action()
      .setImageUri(docker)
      // TODO shouldn't this be using the job shell?
      .setCommands(List("/bin/bash", scriptContainerPath).asJava)
      .setMounts(mounts.asJava)
      .setEntrypoint("")
  }

  def gsutilAsText(command: String*)(mounts: List[Mount] = List.empty, flags: List[ActionFlag] = List.empty): Action = {
    gsutil(List("-h", ContentTypeTextHeader) ++ command.toList: _*)(mounts, flags)
  }

  def gsutil(command: String*)(mounts: List[Mount] = List.empty, flags: List[ActionFlag] = List.empty, description: Option[String] = None): Action = {
    cloudSdkAction
      .setCommands((List("gsutil") ++ command.toList).asJava)
      .withFlags(flags)
      .setMounts(mounts.asJava)
      .setLabels(description.map("description" -> _).toMap.asJava)
  }

  def delocalize(fileOutput: PipelinesApiFileOutput, mounts: List[Mount], projectId: String) = {
    // The command String runs in Bourne shell to get the conditional logic for optional outputs so shell metacharacters in filenames must be escaped.
    val List(containerPath, cloudPath) = List(fileOutput.containerPath.pathAsString, fileOutput.cloudPath) map ESCAPE_XSI.translate

    val copy = s"gsutil -u $projectId cp $containerPath $cloudPath"
    lazy val copyOnlyIfExists = s"if [[ -r $containerPath ]]; then $copy; fi"

    cloudSdkAction
      .withCommand("/bin/sh", "-c", if (fileOutput.optional) copyOnlyIfExists else copy)
      .withFlags(List(ActionFlag.AlwaysRun))
      .withMounts(mounts)
      .withLabels(Map("description" -> "delocalizing"))
  }
}
