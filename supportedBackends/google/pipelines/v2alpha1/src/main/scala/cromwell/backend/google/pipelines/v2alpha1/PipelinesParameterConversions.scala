package cromwell.backend.google.pipelines.v2alpha1

import cats.data.NonEmptyList
import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import com.typesafe.config.ConfigFactory
import cromwell.backend.google.pipelines.common.PipelinesApiAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands._
import cromwell.backend.google.pipelines.v2alpha1.api.{ActionBuilder, ActionFlag}
import cromwell.filesystems.demo.dos.DemoDosPath
import cromwell.filesystems.http.HttpPath
import cromwell.filesystems.sra.SraPath
import simulacrum.typeclass

import scala.language.implicitConversions
@typeclass trait ToParameter[A <: PipelinesParameter] {
  def toActions(p: A, mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration): NonEmptyList[Action]
  def toMount(p: A): Mount = {
    new Mount()
      .setDisk(p.mount.name)
      .setPath(p.mount.mountPoint.pathAsString)
  }
}

trait PipelinesParameterConversions {
  implicit val fileInputToParameter = new ToParameter[PipelinesApiFileInput] {
    override def toActions(fileInput: PipelinesApiFileInput, mounts: List[Mount])
                          (implicit retryPolicy: LocalizationConfiguration): NonEmptyList[Action] = {
      lazy val config = ConfigFactory.load

      val labels = ActionBuilder.parameterLabels(fileInput)
      val describeAction = ActionBuilder.describeParameter(fileInput, labels)
      val localizationAction = fileInput.cloudPath match {
        case _: DemoDosPath =>
          import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands.ShellPath
          import collection.JavaConverters._

          val demoDosDockerImage = config.getString("demo.dos.localization.docker-image")
          val demoDosCommandTemplate = config.getString("demo.dos.localization.command-template")
          val demoDosMarthaUrl = config.getString("demo.dos.martha.url")
          val demoDosCommand = demoDosCommandTemplate
            .replace(s"$${dosPath}", fileInput.cloudPath.escape)
            .replace(s"$${containerPath}", fileInput.containerPath.escape)
          val marthaEnv = Map("MARTHA_URL" -> demoDosMarthaUrl)
          ActionBuilder
            .withImage(demoDosDockerImage)
            .withCommand("/bin/sh", "-c", demoDosCommand)
            .withMounts(mounts)
            .setEnvironment(marthaEnv.asJava)
            .withLabels(labels)
            .setEntrypoint("")
        case sraPath: SraPath =>
          val sraConfig = config.getConfig("filesystems.sra")

          def getString(key: String): Option[String] = {
            if (sraConfig.hasPath(key)) {
              Some(sraConfig.getString(key))
            } else {
              None
            }
          }

          val image = getString("docker-image") getOrElse "fusera/fusera:alpine"
          val (createNgc, ngcArgs) = getString("ngc") match {
            case Some(ngc) => (s"echo $ngc | base64 -d > /tmp/sra.ngc", "-n /tmp/sra.ngc")
            case None => ("", "")
          }
          val mountpoint = s"/cromwell_root/sra-${sraPath.accession}"
          val runFusera = s"fusera mount $ngcArgs -a ${sraPath.accession} $mountpoint"
          ActionBuilder
            .withImage(image)
            .withCommand("/bin/sh", "-c", s"$createNgc; mkdir $mountpoint; $runFusera")
            .withMounts(mounts)
            .withFlags(List(ActionFlag.RunInBackground, ActionFlag.EnableFuse))
        case _: HttpPath =>
          val dockerImage = "google/cloud-sdk:slim"
          val command = s"curl --silent --create-dirs --output ${fileInput.containerPath} ${fileInput.cloudPath}"
          ActionBuilder
            .withImage(dockerImage)
            .withCommand("/bin/sh", "-c", command)
            .withMounts(mounts)
            .withLabels(labels)
            .setEntrypoint("")
        case _ => cloudSdkShellAction(localizeFile(fileInput.cloudPath, fileInput.containerPath))(mounts, labels = labels)
      }

      NonEmptyList.of(describeAction, localizationAction)
    }
  }

  implicit val directoryInputToParameter = new ToParameter[PipelinesApiDirectoryInput] {
    override def toActions(directoryInput: PipelinesApiDirectoryInput, mounts: List[Mount])
                          (implicit retryPolicy: LocalizationConfiguration): NonEmptyList[Action] = {
      val labels = ActionBuilder.parameterLabels(directoryInput)
      val describeAction = ActionBuilder.describeParameter(directoryInput, labels)
      val localizationAction = cloudSdkShellAction(
        localizeDirectory(directoryInput.cloudPath, directoryInput.containerPath)
      )(mounts = mounts, labels = labels)
      NonEmptyList.of(describeAction, localizationAction)
    }
  }

  implicit val fileOutputToParameter = new ToParameter[PipelinesApiFileOutput] {
    override def toActions(fileOutput: PipelinesApiFileOutput, mounts: List[Mount])
                          (implicit retryPolicy: LocalizationConfiguration): NonEmptyList[Action] = {
      // If the output is a "secondary file", it actually could be a directory but we won't know before runtime.
      // The fileOrDirectory method will generate a command that can cover both cases
      val copy = if (fileOutput.secondary)
        delocalizeFileOrDirectory(fileOutput.containerPath, fileOutput.cloudPath, fileOutput.contentType)
      else
        delocalizeFile(fileOutput.containerPath, fileOutput.cloudPath, fileOutput.contentType)

      lazy val copyOnlyIfExists = ifExist(fileOutput.containerPath) { copy }

      val copyCommand = if (fileOutput.optional || fileOutput.secondary) copyOnlyIfExists else copy

      val labels = ActionBuilder.parameterLabels(fileOutput)
      val describeAction = ActionBuilder.describeParameter(fileOutput, labels)

      val delocalizationAction = cloudSdkShellAction(
        copyCommand
      )(mounts = mounts, flags = List(ActionFlag.AlwaysRun), labels = labels)

      // If the file should be uploaded periodically, create 2 actions, a background one with periodic upload, and a normal one
      // that will run at the end and make sure we get the most up to date version of the file
      fileOutput.uploadPeriod match {
        case Some(period) =>
          val periodicLabels = labels collect {
            case (key, _) if key == Key.Tag => key -> Value.Background
            case (key, value) => key -> value
          }
          val periodic = cloudSdkShellAction(
            every(period) { copyCommand }
          )(mounts = mounts, flags = List(ActionFlag.RunInBackground), labels = periodicLabels)

          NonEmptyList.of(describeAction, delocalizationAction, periodic)
        case None => NonEmptyList.of(describeAction, delocalizationAction)
      }
    }
  }

  implicit val directoryOutputToParameter = new ToParameter[PipelinesApiDirectoryOutput] {
    override def toActions(directoryOutput: PipelinesApiDirectoryOutput, mounts: List[Mount])
                          (implicit localizationConfiguration: LocalizationConfiguration): NonEmptyList[Action] = {
      val labels = ActionBuilder.parameterLabels(directoryOutput)
      val describeAction = ActionBuilder.describeParameter(directoryOutput, labels)
      val delocalizationAction = cloudSdkShellAction(
        delocalizeDirectory(directoryOutput.containerPath, directoryOutput.cloudPath, None)
      )(mounts = mounts, flags = List(ActionFlag.AlwaysRun), labels = labels)
      NonEmptyList.of(describeAction, delocalizationAction)
    }
  }

  implicit val inputToParameter = new ToParameter[PipelinesApiInput] {
    override def toActions(p: PipelinesApiInput, mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration) = p match {
      case fileInput: PipelinesApiFileInput => fileInputToParameter.toActions(fileInput, mounts)
      case directoryInput: PipelinesApiDirectoryInput => directoryInputToParameter.toActions(directoryInput, mounts)
    }
  }

  implicit val outputToParameter = new ToParameter[PipelinesApiOutput] {
    override def toActions(p: PipelinesApiOutput, mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration) = p match {
      case fileOutput: PipelinesApiFileOutput => fileOutputToParameter.toActions(fileOutput, mounts)
      case directoryOutput: PipelinesApiDirectoryOutput => directoryOutputToParameter.toActions(directoryOutput, mounts)
    }
  }
}
