package cromwell.backend.google.pipelines.v2alpha1

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionFlag
import simulacrum.typeclass

import scala.language.implicitConversions

@typeclass trait ToParameter[A <: PipelinesParameter] {
  def toActions(p: A, mounts: List[Mount], projectId: String): List[Action]
  def toMount(p: A): Mount = {
    new Mount()
      .setDisk(p.mount.name)
      .setPath(p.mount.mountPoint.pathAsString)
  }
}

trait PipelinesParameterConversions {
  implicit val fileInputToParameter = new ToParameter[PipelinesApiFileInput] {
    override def toActions(fileInput: PipelinesApiFileInput, mounts: List[Mount], projectId: String) = {
      List(gsutil("cp", fileInput.cloudPath.pathAsString, fileInput.containerPath.pathAsString)(mounts, description = Option("localizing")))
    }
  }

  implicit val directoryInputToParameter = new ToParameter[PipelinesApiDirectoryInput] {
    override def toActions(directoryInput: PipelinesApiDirectoryInput, mounts: List[Mount], projectId: String) = {
      // rsync need the target directory to exist already, so create it beforehand
      val mkdirAction = cloudSdkAction
        .withCommand("mkdir", "-p", directoryInput.containerPath.pathAsString)
        .withMounts(mounts)

      val gsutilAction = gsutil("-m", "rsync", "-r", directoryInput.cloudPath.pathAsString, directoryInput.containerPath.pathAsString)(mounts, description = Option("localizing"))

      List(mkdirAction, gsutilAction)
    }
  }

  implicit val fileOutputToParameter = new ToParameter[PipelinesApiFileOutput] {
    override def toActions(fileOutput: PipelinesApiFileOutput, mounts: List[Mount], projectId: String) = {
      List(delocalizeFile(fileOutput, mounts, projectId))
    }
  }

  implicit val directoryOutputToParameter = new ToParameter[PipelinesApiDirectoryOutput] {
    override def toActions(directoryOutput: PipelinesApiDirectoryOutput, mounts: List[Mount], projectId: String) = {
      List(gsutil("-m", "rsync", "-r", directoryOutput.containerPath.pathAsString, directoryOutput.cloudPath.pathAsString)(mounts, List(ActionFlag.AlwaysRun), description = Option("delocalizing")))
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
