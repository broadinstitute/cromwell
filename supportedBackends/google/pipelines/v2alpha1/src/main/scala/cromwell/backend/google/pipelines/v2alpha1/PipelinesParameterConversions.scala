package cromwell.backend.google.pipelines.v2alpha1

import cats.data.NonEmptyList
import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import com.typesafe.config.ConfigFactory
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands._
import cromwell.backend.google.pipelines.v2alpha1.api.{ActionBuilder, ActionFlag}
import simulacrum.typeclass

import scala.language.implicitConversions
@typeclass trait ToParameter[A <: PipelinesParameter] {
  def toActions(p: A, mounts: List[Mount]): NonEmptyList[Action]
  def toMount(p: A): Mount = {
    new Mount()
      .setDisk(p.mount.name)
      .setPath(p.mount.mountPoint.pathAsString)
  }
}

trait PipelinesParameterConversions {
  implicit val fileInputToParameter = new ToParameter[PipelinesApiFileInput] {
    override def toActions(fileInput: PipelinesApiFileInput, mounts: List[Mount]) = NonEmptyList.of {
      val labels = Map(
        Key.Tag -> Value.Localization,
        Key.InputName -> fileInput.name
      )

      if (fileInput.cloudPath.pathAsString.startsWith("dos://")) {
        import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands.ShellPath
        val config = ConfigFactory.load
        val demoDosDockerImage = config.getString("demo.dos.docker-image")
        val demoDosCommandTemplate = config.getString("demo.dos.command-template")
        val demoDosCommand = demoDosCommandTemplate
          .replace(s"$${cloudPath}", fileInput.cloudPath.escape)
          .replace(s"$${containerPath}", fileInput.containerPath.escape)
        ActionBuilder
          .withImage(demoDosDockerImage)
          .withCommand("/bin/sh", "-c", demoDosCommand)
          .withMounts(mounts)
          .withLabels(labels)
          .setEntrypoint("")
      } else {
        cloudSdkShellAction(localizeFile(fileInput.cloudPath, fileInput.containerPath))(mounts, labels = labels)
      }

    }
  }

  implicit val directoryInputToParameter = new ToParameter[PipelinesApiDirectoryInput] {
    override def toActions(directoryInput: PipelinesApiDirectoryInput, mounts: List[Mount]) =  NonEmptyList.of {
      cloudSdkShellAction(
        localizeDirectory(directoryInput.cloudPath, directoryInput.containerPath)
      )(mounts, labels =  Map(Key.Tag -> Value.Localization))
    }
  }

  implicit val fileOutputToParameter = new ToParameter[PipelinesApiFileOutput] {
    override def toActions(fileOutput: PipelinesApiFileOutput, mounts: List[Mount]) = {
      // If the output is a "secondary file", it actually could be a directory but we won't know before runtime.
      // The fileOrDirectory method will generate a command that can cover both cases
      val copy = if (fileOutput.secondary)
        delocalizeFileOrDirectory(fileOutput.containerPath, fileOutput.cloudPath, fileOutput.contentType)
      else
        delocalizeFile(fileOutput.containerPath, fileOutput.cloudPath, fileOutput.contentType)

      lazy val copyOnlyIfExists = ifExist(fileOutput.containerPath) { copy }

      val copyCommand = if (fileOutput.optional || fileOutput.secondary) copyOnlyIfExists else copy
      
      val delocalizationAction = cloudSdkShellAction(
        copyCommand
      )(mounts = mounts, flags = List(ActionFlag.AlwaysRun), labels = Map(Key.Tag -> Value.Delocalization))

      // If the file should be uploaded periodically, create 2 actions, a background one with periodic upload, and a normal one
      // that will run at the end and make sure we get the most up to date version of the file
      fileOutput.uploadPeriod match {
        case Some(period) =>
          val periodic = cloudSdkShellAction(
            every(period) { copyCommand }
          )(mounts = mounts, flags = List(ActionFlag.RunInBackground), labels = Map(Key.Tag -> Value.Background))

          NonEmptyList.of(delocalizationAction, periodic)
        case None => NonEmptyList.of(delocalizationAction)
      }
    }
  }

  implicit val directoryOutputToParameter = new ToParameter[PipelinesApiDirectoryOutput] {
    override def toActions(directoryOutput: PipelinesApiDirectoryOutput, mounts: List[Mount]) = NonEmptyList.of {
      cloudSdkShellAction(
        delocalizeDirectory(directoryOutput.containerPath, directoryOutput.cloudPath, None)
      )(mounts, List(ActionFlag.AlwaysRun), labels =  Map(Key.Tag -> Value.Delocalization))
    }
  }

  implicit val inputToParameter = new ToParameter[PipelinesApiInput] {
    override def toActions(p: PipelinesApiInput, mounts: List[Mount]) = p match {
      case fileInput: PipelinesApiFileInput => fileInputToParameter.toActions(fileInput, mounts)
      case directoryInput: PipelinesApiDirectoryInput => directoryInputToParameter.toActions(directoryInput, mounts)
    }
  }

  implicit val outputToParameter = new ToParameter[PipelinesApiOutput] {
    override def toActions(p: PipelinesApiOutput, mounts: List[Mount]) = p match {
      case fileOutput: PipelinesApiFileOutput => fileOutputToParameter.toActions(fileOutput, mounts)
      case directoryOutput: PipelinesApiDirectoryOutput => directoryOutputToParameter.toActions(directoryOutput, mounts)
    }
  }
}
