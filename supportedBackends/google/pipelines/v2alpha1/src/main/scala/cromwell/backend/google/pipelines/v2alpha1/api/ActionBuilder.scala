package cromwell.backend.google.pipelines.v2alpha1.api

import akka.http.scaladsl.model.ContentTypes
import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Gsutil.ContentTypeTextHeader
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionFlag.ActionFlag
import mouse.all._

import scala.collection.JavaConverters._

/**
  * Utility singleton to create high level actions.
  */
object ActionBuilder {
  object Labels {
    object Key {
      /**
        * Very short description of the action
        */
      val Tag = "tag"
      val InputName = "inputName"
    }
    object Value {
      val UserAction = "UserAction"
      val Localization = "Localization"
      val Delocalization = "Delocalization"
    }
  }

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

  def userAction(docker: String, scriptContainerPath: String, mounts: List[Mount], jobShell: String): Action = {
    new Action()
      .setImageUri(docker)
      .setCommands(List(jobShell, scriptContainerPath).asJava)
      .setMounts(mounts.asJava)
      .setEntrypoint("")
      .setLabels(Map(Key.Tag -> Value.UserAction).asJava)
  }

  def gsutilAsText(command: String*)(mounts: List[Mount] = List.empty, flags: List[ActionFlag] = List.empty, labels: Map[String, String] = Map.empty): Action = {
    gsutil(List("-h", ContentTypeTextHeader) ++ command.toList: _*)(mounts, flags)
  }

  def gsutil(command: String*)(mounts: List[Mount] = List.empty, flags: List[ActionFlag] = List.empty, labels: Map[String, String] = Map.empty): Action = {
    cloudSdkAction
      .setCommands((List("gsutil") ++ command.toList).asJava)
      .withFlags(flags)
      .setMounts(mounts.asJava)
      .setLabels(labels.asJava)
  }

  def cloudSdkBashAction(shellCommand: String)(mounts: List[Mount] = List.empty, flags: List[ActionFlag] = List.empty, labels: Map[String, String] = Map.empty): Action = 
    cloudSdkAction
      .withCommand("/bin/sh", "-c", shellCommand)
      .withFlags(flags)
      .withMounts(mounts)
      .withLabels(labels)
}
