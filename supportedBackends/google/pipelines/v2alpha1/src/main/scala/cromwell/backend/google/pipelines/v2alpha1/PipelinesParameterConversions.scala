package cromwell.backend.google.pipelines.v2alpha1

import cats.data.NonEmptyList
import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionFlag
import simulacrum.typeclass

import scala.language.implicitConversions
@typeclass trait ToParameter[A <: PipelinesParameter] {
  def toActions(p: A, mounts: List[Mount], projectId: String): NonEmptyList[Action]
  def toMount(p: A): Mount = {
    new Mount()
      .setDisk(p.mount.name)
      .setPath(p.mount.mountPoint.pathAsString)
  }
}

trait PipelinesParameterConversions {
  implicit val fileInputToParameter = new ToParameter[PipelinesApiFileInput] {
    override def toActions(fileInput: PipelinesApiFileInput, mounts: List[Mount], projectId: String) = NonEmptyList.of {
      val labels = Map(
        Key.Tag -> Value.Localization,
        Key.InputName -> fileInput.name
      )
      cloudSdkBashAction(localizeFile(fileInput.cloudPath, fileInput.containerPath))(mounts, labels = labels)
    }
  }

  implicit val directoryInputToParameter = new ToParameter[PipelinesApiDirectoryInput] {
    override def toActions(directoryInput: PipelinesApiDirectoryInput, mounts: List[Mount], projectId: String) =  NonEmptyList.of {
      cloudSdkBashAction(
        localizeDirectory(directoryInput.cloudPath, directoryInput.containerPath)
      )(mounts, labels =  Map(Key.Tag -> Value.Localization))
    }
  }

  implicit val fileOutputToParameter = new ToParameter[PipelinesApiFileOutput] {
    override def toActions(fileOutput: PipelinesApiFileOutput, mounts: List[Mount], projectId: String) = NonEmptyList.of {
      // If the output is a "secondary file", it actually could be a directory but we won't know before runtime.
      // The fileOrDirectory method will generate a command that can cover both cases
      val copy = if (fileOutput.secondary)
        delocalizeFileOrDirectory(fileOutput.containerPath, fileOutput.cloudPath)
      else
        delocalizeFile(fileOutput.containerPath, fileOutput.cloudPath)

      lazy val copyOnlyIfExists = ifExist(fileOutput.containerPath) { copy }

      cloudSdkBashAction(
        if (fileOutput.optional || fileOutput.secondary) copyOnlyIfExists else copy
      )(mounts = mounts, flags = List(ActionFlag.AlwaysRun), labels = Map(Key.Tag -> Value.Delocalization))
    }
  }

  implicit val directoryOutputToParameter = new ToParameter[PipelinesApiDirectoryOutput] {
    override def toActions(directoryOutput: PipelinesApiDirectoryOutput, mounts: List[Mount], projectId: String) = NonEmptyList.of {
      cloudSdkBashAction(
        delocalizeDirectory(directoryOutput.containerPath, directoryOutput.cloudPath)
      )(mounts, List(ActionFlag.AlwaysRun), labels =  Map(Key.Tag -> Value.Delocalization))
    }
  }

  implicit val inputToParameter = new ToParameter[PipelinesApiInput] {
    override def toActions(p: PipelinesApiInput, mounts: List[Mount], projectId: String) = p match {
      case fileInput: PipelinesApiFileInput => fileInputToParameter.toActions(fileInput, mounts, projectId)
      case directoryInput: PipelinesApiDirectoryInput => directoryInputToParameter.toActions(directoryInput, mounts, projectId)
    }
  }

  implicit val outputToParameter = new ToParameter[PipelinesApiOutput] {
    override def toActions(p: PipelinesApiOutput, mounts: List[Mount], projectId: String) = p match {
      case fileOutput: PipelinesApiFileOutput => fileOutputToParameter.toActions(fileOutput, mounts, projectId)
      case directoryOutput: PipelinesApiDirectoryOutput => directoryOutputToParameter.toActions(directoryOutput, mounts, projectId)
    }
  }
}
