package cromwell.backend.google.pipelines.batch.runnable

import cloud.nio.impl.drs.DrsConfig
import com.typesafe.config.ConfigFactory
import com.google.cloud.batch.v1.{Runnable, Volume}

import cromwell.backend.google.pipelines.batch.runnable.RunnableCommands.localizeFile
import cromwell.backend.google.pipelines.common.action.ActionLabels._
import cromwell.backend.google.pipelines.batch.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.batch.GcpBatchFileInput
import cromwell.backend.google.pipelines.batch.GcpBatchJobPaths._
import cromwell.backend.google.pipelines.batch.GcpBatchRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.batch.BatchParameterConversions._
import cromwell.backend.google.pipelines.batch.runnable.RunnableBuilder.{EnhancedRunnable, cloudSdkShellAction}
import cromwell.backend.google.pipelines.batch.runnable.RunnableCommands._


import cromwell.core.path.Path
import cromwell.filesystems.drs.DrsPath
import scala.jdk.CollectionConverters._


trait Localization {


  def localizeActions(createPipelineParameters: CreatePipelineParameters, volumes: List[Volume])
                     (implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Runnable] = {
    val localizationLabel = Map(Key.Tag -> Value.Localization)

    val gcsTransferLibraryContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsTransferLibraryName)
    val localizeGcsTransferLibrary = cloudSdkShellAction(localizeFile(
      cloudPath = createPipelineParameters.cloudCallRoot / GcsTransferLibraryName,
      containerPath = gcsTransferLibraryContainerPath))(volumes = volumes, labels = localizationLabel)

    val gcsLocalizationContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsLocalizationScriptName)
    val localizeGcsLocalizationScript = cloudSdkShellAction(localizeFile(
      cloudPath = createPipelineParameters.cloudCallRoot / GcsLocalizationScriptName,
      containerPath = gcsLocalizationContainerPath))(volumes = volumes, labels = localizationLabel)

    val gcsDelocalizationContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsDelocalizationScriptName)
    val localizeGcsDelocalizationScript = cloudSdkShellAction(localizeFile(
      cloudPath = createPipelineParameters.cloudCallRoot / GcsDelocalizationScriptName,
      containerPath = gcsDelocalizationContainerPath))(volumes = volumes, labels = localizationLabel)

    val runGcsLocalizationScript = cloudSdkShellAction(
      s"/bin/bash $gcsLocalizationContainerPath")(volumes = volumes, labels = localizationLabel)

    val drsInputs: List[DrsPath] = createPipelineParameters.inputOutputParameters.fileInputParameters.collect {
      case GcpBatchFileInput(_, drsPath: DrsPath, _, _) => drsPath
    }

    val drsLocalizationActions = if (drsInputs.nonEmpty) {
      val drsLocalizationManifestContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(DrsLocalizationManifestName)
      val localizeDrsLocalizationManifest = cloudSdkShellAction(localizeFile(
        cloudPath = createPipelineParameters.cloudCallRoot / DrsLocalizationManifestName,
        containerPath = drsLocalizationManifestContainerPath))(volumes = volumes, labels = localizationLabel)

      // Requester pays project id is stored on each DrsPath, but will be the same for all DRS inputs to a
      // particular workflow because it's determined by the Google project set in workflow options.
      val requesterPaysProjectId: Option[String] = drsInputs.flatMap(_.requesterPaysProjectIdOption).headOption
      val runDrsLocalization = Localization.drsAction(drsLocalizationManifestContainerPath, volumes, localizationLabel, requesterPaysProjectId)
      List(localizeDrsLocalizationManifest, runDrsLocalization)
    } else List[Runnable]()

    // Any "classic" PAPI v2 one-at-a-time localizations for non-GCS inputs.
    val singletonLocalizations = createPipelineParameters.inputOutputParameters.fileInputParameters.flatMap(_.toRunnables(volumes).toList)

    val localizations =
      localizeGcsTransferLibrary ::
        localizeGcsLocalizationScript :: runGcsLocalizationScript ::
        drsLocalizationActions :::
        localizeGcsDelocalizationScript ::
        singletonLocalizations

    RunnableBuilder.annotateTimestampedActions("localization", Value.Localization)(localizations)
  }

}
object Localization {

  def drsAction(manifestPath: Path,
                volumes: List[Volume],
                labels: Map[String, String],
                requesterPaysProjectId: Option[String]
               ): Runnable = {
    val config = ConfigFactory.load
    val marthaConfig = config.getConfig("filesystems.drs.global.config.martha")
    val drsConfig = DrsConfig.fromConfig(marthaConfig)
    val drsDockerImage = config.getString("drs.localization.docker-image")

    val manifestArg = List("-m", manifestPath.pathAsString)
    val requesterPaysArg = requesterPaysProjectId.map(r => List("-r", r)).getOrElse(List.empty)
    val drsCommand = manifestArg ++ requesterPaysArg

    val marthaEnv = DrsConfig.toEnv(drsConfig)
    RunnableBuilder
      .withImage(drsDockerImage)
      .withCommand(drsCommand: _*)
      .withMounts(volumes)
      .setEnvironment(marthaEnv.asJava)
      .withLabels(labels)
  }
}
