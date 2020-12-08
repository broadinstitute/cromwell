package cromwell.backend.google.pipelines.v2beta

import cloud.nio.impl.drs.{DrsCloudNioFileSystemProvider, DrsConfig}
import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import com.typesafe.config.ConfigFactory
import cromwell.backend.google.pipelines.common.action.ActionLabels._
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.v2beta.api.ActionBuilder
import cromwell.backend.google.pipelines.v2beta.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2beta.api.ActionCommands._
import cromwell.filesystems.drs.DrsPath
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.http.HttpPath
import cromwell.filesystems.sra.SraPath
import simulacrum.typeclass

@typeclass trait ToParameter[A <: PipelinesParameter] {
  def toActions(p: A, mounts: List[Mount])(implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Action]
}

trait PipelinesParameterConversions {
  implicit val fileInputToParameter: ToParameter[PipelinesApiFileInput] = new ToParameter[PipelinesApiFileInput] {
    override def toActions(fileInput: PipelinesApiFileInput, mounts: List[Mount])
                          (implicit retryPolicy: GcsTransferConfiguration): List[Action] = {
      lazy val config = ConfigFactory.load

      val labels = ActionBuilder.parameterLabels(fileInput)
      fileInput.cloudPath match {
        case drsPath: DrsPath =>

          import collection.JavaConverters._

          val drsFileSystemProvider = drsPath.drsPath.getFileSystem.provider.asInstanceOf[DrsCloudNioFileSystemProvider]

          val drsDockerImage = config.getString("drs.localization.docker-image")
          // Note: Don't ShellPath.escape the paths as we are directly invoking the localizer and NOT launching a shell.
          val drsCommand =
            List(fileInput.cloudPath.pathAsString, fileInput.containerPath.pathAsString) ++
              drsPath.requesterPaysProjectIdOption.toList
          val marthaEnv = DrsConfig.toEnv(drsFileSystemProvider.drsConfig)
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
            .withEntrypointCommand("/bin/sh", "-c", s"$createNgc; mkdir $mountpoint; $runFusera")
            .withMounts(mounts)
            .withRunInBackground(true)
            .withEnableFuse(true)
          List(ActionBuilder.describeParameter(fileInput, labels), localizationAction)
        case _: HttpPath =>
          val command = s"curl --silent --create-dirs --output ${fileInput.containerPath} ${fileInput.cloudPath}"
          val localizationAction = ActionBuilder.cloudSdkShellAction(command)(mounts = mounts, labels = labels)
          List(ActionBuilder.describeParameter(fileInput, labels), localizationAction)
        case _: GcsPath =>
          // GCS paths will be localized with a separate localization script.
          Nil
      }
    }
  }

  implicit val directoryInputToParameter: ToParameter[PipelinesApiDirectoryInput] =
    new ToParameter[PipelinesApiDirectoryInput] {
    override def toActions(directoryInput: PipelinesApiDirectoryInput, mounts: List[Mount])
                          (implicit retryPolicy: GcsTransferConfiguration): List[Action] = {
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

  implicit val fileOutputToParameter: ToParameter[PipelinesApiFileOutput] = new ToParameter[PipelinesApiFileOutput] {
    override def toActions(fileOutput: PipelinesApiFileOutput, mounts: List[Mount])
                          (implicit retryPolicy: GcsTransferConfiguration): List[Action] = {

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

          val delocalizationAction = cloudSdkShellAction(copyCommand)(mounts = mounts, labels = labels).withAlwaysRun(true)
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
          )(mounts = mounts, labels = periodicLabels).withRunInBackground(true)

          finalDelocalizationActions :+ periodic
        case None => finalDelocalizationActions
      }
    }
  }

  implicit val directoryOutputToParameter: ToParameter[PipelinesApiDirectoryOutput] =
    new ToParameter[PipelinesApiDirectoryOutput] {
    override def toActions(directoryOutput: PipelinesApiDirectoryOutput, mounts: List[Mount])
                          (implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Action] = {
      directoryOutput.cloudPath match {
        case _: GcsPath => Nil // GCS paths will be delocalized with a separate delocalization script.
        case _ =>
          val labels = ActionBuilder.parameterLabels(directoryOutput)
          val describeAction = ActionBuilder.describeParameter(directoryOutput, labels)
          val delocalizationAction = cloudSdkShellAction(
            delocalizeDirectory(directoryOutput.containerPath, directoryOutput.cloudPath, None)
          )(mounts = mounts, labels = labels).withAlwaysRun(true)
          List(describeAction, delocalizationAction)
      }
    }
  }

  implicit val inputToParameter: ToParameter[PipelinesApiInput] = new ToParameter[PipelinesApiInput] {
    override def toActions(p: PipelinesApiInput, mounts: List[Mount])
                          (implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Action] = p match {
      case fileInput: PipelinesApiFileInput => fileInputToParameter.toActions(fileInput, mounts)
      case directoryInput: PipelinesApiDirectoryInput => directoryInputToParameter.toActions(directoryInput, mounts)
    }
  }

  implicit val outputToParameter: ToParameter[PipelinesApiOutput] = new ToParameter[PipelinesApiOutput] {
    override def toActions(p: PipelinesApiOutput, mounts: List[Mount])
                          (implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Action] = p match {
      case fileOutput: PipelinesApiFileOutput => fileOutputToParameter.toActions(fileOutput, mounts)
      case directoryOutput: PipelinesApiDirectoryOutput => directoryOutputToParameter.toActions(directoryOutput, mounts)
    }
  }
}
