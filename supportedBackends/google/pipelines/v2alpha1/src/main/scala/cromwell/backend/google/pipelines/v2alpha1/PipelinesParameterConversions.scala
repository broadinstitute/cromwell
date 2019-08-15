package cromwell.backend.google.pipelines.v2alpha1

import cloud.nio.impl.drs.DrsCloudNioFileSystemProvider
import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import com.typesafe.config.ConfigFactory
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands._
import cromwell.backend.google.pipelines.v2alpha1.api.{ActionBuilder, ActionFlag}
import cromwell.filesystems.drs.DrsPath
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.http.HttpPath
import cromwell.filesystems.sra.SraPath
import simulacrum.typeclass

import scala.language.implicitConversions
@typeclass trait ToParameter[A <: PipelinesParameter] {
  def toActions(p: A, mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration): List[Action]
}

trait PipelinesParameterConversions {
  implicit val fileInputToParameter = new ToParameter[PipelinesApiFileInput] {
    override def toActions(fileInput: PipelinesApiFileInput, mounts: List[Mount])
                          (implicit retryPolicy: LocalizationConfiguration): List[Action] = {
      lazy val config = ConfigFactory.load

      val labels = ActionBuilder.parameterLabels(fileInput)
      fileInput.cloudPath match {
        case drsPath: DrsPath =>
          import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands.ShellPath
          import collection.JavaConverters._

          val drsFileSystemProvider = drsPath.drsPath.getFileSystem.provider.asInstanceOf[DrsCloudNioFileSystemProvider]

          val drsDockerImage = config.getString("drs.localization.docker-image")
          val drsMarthaUrl = drsFileSystemProvider.config.getString("martha.url")
          val drsCommand = List(fileInput.cloudPath.escape, fileInput.containerPath.escape) ++ drsPath.requesterPaysProjectIdOption.toList
          val marthaEnv = Map("MARTHA_URL" -> drsMarthaUrl)
          val localizationAction = ActionBuilder
            .withImage(drsDockerImage)
            .withCommand(drsCommand: _*)
            .withMounts(mounts)
            .setEnvironment(marthaEnv.asJava)
            .withLabels(labels)
          List(ActionBuilder.describeParameter(fileInput, labels), localizationAction)
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
          val localizationAction = ActionBuilder
            .withImage(image)
            .withCommand("/bin/sh", "-c", s"$createNgc; mkdir $mountpoint; $runFusera")
            .withMounts(mounts)
            .withFlags(List(ActionFlag.RunInBackground, ActionFlag.EnableFuse))
          List(ActionBuilder.describeParameter(fileInput, labels), localizationAction)
        case _: HttpPath =>
          val dockerImage = "google/cloud-sdk:slim"
          val command = s"curl --silent --create-dirs --output ${fileInput.containerPath} ${fileInput.cloudPath}"
          val localizationAction = ActionBuilder
            .withImage(dockerImage)
            .withCommand("/bin/sh", "-c", command)
            .withMounts(mounts)
            .withLabels(labels)
            .setEntrypoint("")
          List(ActionBuilder.describeParameter(fileInput, labels), localizationAction)
        case _: GcsPath =>
          // GCS paths will be localized with a separate localization script.
          Nil
      }
    }
  }

  implicit val directoryInputToParameter = new ToParameter[PipelinesApiDirectoryInput] {
    override def toActions(directoryInput: PipelinesApiDirectoryInput, mounts: List[Mount])
                          (implicit retryPolicy: LocalizationConfiguration): List[Action] = {
      directoryInput.cloudPath match {
        case _: GcsPath => Nil // GCS paths will be localized with a separate localization script.
        case _ =>
          val labels = ActionBuilder.parameterLabels(directoryInput)
          val describeAction = ActionBuilder.describeParameter(directoryInput, labels)
          val localizationAction = cloudSdkShellAction(
            localizeDirectory(directoryInput.cloudPath, directoryInput.containerPath)
          )(mounts = mounts, labels = labels)
          List(describeAction, localizationAction)
      }
    }
  }

  implicit val fileOutputToParameter = new ToParameter[PipelinesApiFileOutput] {
    override def toActions(fileOutput: PipelinesApiFileOutput, mounts: List[Mount])
                          (implicit retryPolicy: LocalizationConfiguration): List[Action] = {

      // If the output is a "secondary file", it actually could be a directory but we won't know before runtime.
      // The fileOrDirectory method will generate a command that can cover both cases
      lazy val copy = if (fileOutput.secondary)
        delocalizeFileOrDirectory(fileOutput.containerPath, fileOutput.cloudPath, fileOutput.contentType)
      else
        delocalizeFile(fileOutput.containerPath, fileOutput.cloudPath, fileOutput.contentType)

      lazy val copyOnlyIfExists = ifExist(fileOutput.containerPath) {
        copy
      }

      lazy val copyCommand = if (fileOutput.optional || fileOutput.secondary) copyOnlyIfExists else copy
      lazy val labels = ActionBuilder.parameterLabels(fileOutput)

      // The delocalization actions to take once the user command has terminated (i.e., the non-periodic uploads).
      val finalDelocalizationActions = fileOutput.cloudPath match {
        case _: GcsPath => Nil // GCS files are delocalized with a separate delocalization script.
        case _ =>
          val describeAction = ActionBuilder.describeParameter(fileOutput, labels)

          val delocalizationAction = cloudSdkShellAction(
            copyCommand
          )(mounts = mounts, flags = List(ActionFlag.AlwaysRun), labels = labels)
          List(describeAction, delocalizationAction)
      }

      fileOutput.uploadPeriod match {
        // If the file should be uploaded periodically, create a background upload action in addition to any normal ones
        // that run at the end to make sure we get the most up to date version of the file.
        case Some(period) =>
          val periodicLabels = labels collect {
            case (key, _) if key == Key.Tag => key -> Value.Background
            case (key, value) => key -> value
          }
          val periodic = cloudSdkShellAction(
            every(period) {
              copyCommand
            }
          )(mounts = mounts, flags = List(ActionFlag.RunInBackground), labels = periodicLabels)

          finalDelocalizationActions :+ periodic
        case None => finalDelocalizationActions
      }
    }
  }

  implicit val directoryOutputToParameter = new ToParameter[PipelinesApiDirectoryOutput] {
    override def toActions(directoryOutput: PipelinesApiDirectoryOutput, mounts: List[Mount])
                          (implicit localizationConfiguration: LocalizationConfiguration): List[Action] = {
      directoryOutput.cloudPath match {
        case _: GcsPath => Nil // GCS paths will be delocalized with a separate delocalization script.
        case _ =>
          val labels = ActionBuilder.parameterLabels(directoryOutput)
          val describeAction = ActionBuilder.describeParameter(directoryOutput, labels)
          val delocalizationAction = cloudSdkShellAction(
            delocalizeDirectory(directoryOutput.containerPath, directoryOutput.cloudPath, None)
          )(mounts = mounts, flags = List(ActionFlag.AlwaysRun), labels = labels)
          List(describeAction, delocalizationAction)
      }
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
