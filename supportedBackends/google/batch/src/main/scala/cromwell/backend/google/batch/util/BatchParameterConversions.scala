package cromwell.backend.google.batch.util

import com.google.cloud.batch.v1.{Runnable, Volume}
//import com.typesafe.config.ConfigFactory
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.batch.models._
import cromwell.backend.google.batch.runnable._
import cromwell.filesystems.drs.DrsPath
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.http.HttpPath
//import cromwell.filesystems.sra.SraPath
import simulacrum.typeclass

@typeclass trait ToParameter[A <: BatchParameter] {
  def toRunnables(p: A, volumes: List[Volume])(implicit
    gcsTransferConfiguration: GcsTransferConfiguration
  ): List[Runnable.Builder]
}

trait GcpBatchParameterConversions {
  import RunnableBuilder._
  import RunnableCommands._
  import RunnableLabels._

  implicit val fileInputToParameter: ToParameter[GcpBatchFileInput] = new ToParameter[GcpBatchFileInput] {
    override def toRunnables(fileInput: GcpBatchFileInput, volumes: List[Volume])(implicit
      retryPolicy: GcsTransferConfiguration
    ): List[Runnable.Builder] = {

      val labels = RunnableBuilder.parameterLabels(fileInput)
      fileInput.cloudPath match {
        case _: HttpPath =>
          val command = s"curl --silent --create-dirs --output ${fileInput.containerPath} ${fileInput.cloudPath}"
          val localizationRunnables =
            RunnableBuilder.cloudSdkShellRunnable(command)(volumes = volumes, labels = labels, flags = List.empty)
          List(RunnableBuilder.describeParameter(fileInput, volumes, labels), localizationRunnables)

        case _: GcsPath =>
          // GCS paths will be localized with a separate localization script.
          Nil
        case _: DrsPath =>
          // DRS paths will be localized with a single call to cromwell-drs-localizer with a manifest
          Nil
      }
    }
  }

  implicit val directoryInputToParameter: ToParameter[GcpBatchDirectoryInput] =
    new ToParameter[GcpBatchDirectoryInput] {
      override def toRunnables(directoryInput: GcpBatchDirectoryInput, volumes: List[Volume])(implicit
        retryPolicy: GcsTransferConfiguration
      ): List[Runnable.Builder] =
        directoryInput.cloudPath match {
          case _: GcsPath => Nil // GCS paths will be localized with a separate localization script.
          case _ =>
            val labels = RunnableBuilder.parameterLabels(directoryInput)
            val describeRunnables = RunnableBuilder.describeParameter(directoryInput, volumes, labels)
            val localizationRunnables = RunnableBuilder.cloudSdkShellRunnable(
              RunnableCommands.localizeDirectory(directoryInput.cloudPath, directoryInput.containerPath)
            )(volumes = volumes, labels = labels, flags = List.empty)
            List(describeRunnables, localizationRunnables)
        }
    }

  implicit val fileOutputToParameter: ToParameter[GcpBatchFileOutput] = new ToParameter[GcpBatchFileOutput] {
    override def toRunnables(fileOutput: GcpBatchFileOutput, volumes: List[Volume])(implicit
      retryPolicy: GcsTransferConfiguration
    ): List[Runnable.Builder] = {

      // If the output is a "secondary file", it actually could be a directory but we won't know before runtime.
      // The fileOrDirectory method will generate a command that can cover both cases
      lazy val copy =
        if (fileOutput.secondary)
          RunnableCommands.delocalizeFileOrDirectory(fileOutput.containerPath,
                                                     fileOutput.cloudPath,
                                                     fileOutput.contentType
          )
        else
          RunnableCommands.delocalizeFile(fileOutput.containerPath, fileOutput.cloudPath, fileOutput.contentType)

      lazy val copyOnlyIfExists = RunnableCommands.ifExist(fileOutput.containerPath) {
        copy
      }

      lazy val copyCommand = if (fileOutput.optional || fileOutput.secondary) copyOnlyIfExists else copy
      lazy val labels = RunnableBuilder.parameterLabels(fileOutput)

      // The delocalization runnables to take once the user command has terminated (i.e., the non-periodic uploads).
      val finalDelocalizationRunnables = fileOutput.cloudPath match {
        case _: GcsPath => Nil // GCS files are delocalized with a separate delocalization script.
        case _ =>
          val describeRunnable = RunnableBuilder.describeParameter(fileOutput, volumes, labels)
          val delocalizationRunnable = RunnableBuilder
            .cloudSdkShellRunnable(copyCommand)(volumes = volumes, labels = labels, flags = List.empty)
            .withAlwaysRun(true)

          List(describeRunnable, delocalizationRunnable)
      }

      fileOutput.uploadPeriod match {
        // If the file should be uploaded periodically, create a background upload runnable in addition to any normal ones
        // that run at the end to make sure we get the most up to date version of the file.
        case Some(period) =>
          val periodicLabels = labels collect {
            case (key, _) if key == Key.Tag => key -> Value.Background
            case (key, value) => key -> value
          }
          val periodic = RunnableBuilder
            .cloudSdkShellRunnable(
              every(period) {
                copyCommand
              }
            )(volumes = volumes, labels = periodicLabels, flags = List.empty)
            .withRunInBackground(true)

          finalDelocalizationRunnables :+ periodic

        case None => finalDelocalizationRunnables
      }
    }
  }

  implicit val directoryOutputToParameter: ToParameter[GcpBatchDirectoryOutput] =
    new ToParameter[GcpBatchDirectoryOutput] {
      override def toRunnables(directoryOutput: GcpBatchDirectoryOutput, volumes: List[Volume])(implicit
        gcsTransferConfiguration: GcsTransferConfiguration
      ): List[Runnable.Builder] =
        directoryOutput.cloudPath match {
          case _: GcsPath => Nil // GCS paths will be delocalized with a separate delocalization script.
          case _ =>
            val labels = RunnableBuilder.parameterLabels(directoryOutput)
            val describeRunnable = RunnableBuilder.describeParameter(directoryOutput, volumes, labels)
            val delocalizationRunnable = RunnableBuilder
              .cloudSdkShellRunnable(
                delocalizeDirectory(directoryOutput.containerPath, directoryOutput.cloudPath, None)
              )(volumes = volumes, labels = labels, flags = List.empty)
              .withAlwaysRun(true)

            List(describeRunnable, delocalizationRunnable)
        }
    }

  implicit val inputToParameter: ToParameter[GcpBatchInput] = new ToParameter[GcpBatchInput] {
    override def toRunnables(p: GcpBatchInput, volumes: List[Volume])(implicit
      gcsTransferConfiguration: GcsTransferConfiguration
    ): List[Runnable.Builder] = p match {
      case fileInput: GcpBatchFileInput => fileInputToParameter.toRunnables(fileInput, volumes)
      case directoryInput: GcpBatchDirectoryInput => directoryInputToParameter.toRunnables(directoryInput, volumes)
    }
  }

  implicit val outputToParameter: ToParameter[GcpBatchOutput] = new ToParameter[GcpBatchOutput] {
    override def toRunnables(p: GcpBatchOutput, volumes: List[Volume])(implicit
      gcsTransferConfiguration: GcsTransferConfiguration
    ): List[Runnable.Builder] = p match {
      case fileOutput: GcpBatchFileOutput => fileOutputToParameter.toRunnables(fileOutput, volumes)
      case directoryOutput: GcpBatchDirectoryOutput => directoryOutputToParameter.toRunnables(directoryOutput, volumes)
    }
  }
}

object GcpBatchParameterConversions extends GcpBatchParameterConversions
