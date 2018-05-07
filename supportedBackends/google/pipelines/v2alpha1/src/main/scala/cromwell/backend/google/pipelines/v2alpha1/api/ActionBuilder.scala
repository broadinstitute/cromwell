package cromwell.backend.google.pipelines.v2alpha1.api

import akka.http.scaladsl.model.ContentTypes
import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Gsutil.ContentTypeTextHeader
import cromwell.backend.google.pipelines.v2alpha1.api.ActionFlag.ActionFlag
import mouse.all._

import scala.collection.JavaConverters._

/**
  * Utility singleton to create high level actions.
  */
object ActionBuilder {
  object Gsutil {
    private val contentTypeText = ContentTypes.`text/plain(UTF-8)`.toString()
    val ContentTypeTextHeader = s"Content-Type: $contentTypeText"
  }

  private val cloudSdkImage = "google/cloud-sdk:alpine"
  def cloudSdkAction: Action = new Action().setImageUri(cloudSdkImage)

  private def javaFlags(flags: List[ActionFlag]) = flags.map(_.toString).asJava
  
  def userAction(docker: String, scriptContainerPath: String, mounts: List[Mount]): Action = {
    new Action()
      .setImageUri(docker)
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
      .setFlags(flags |> javaFlags)
      .setMounts(mounts.asJava)
      .setLabels(description.map("description" -> _).toMap.asJava)
  }
}
