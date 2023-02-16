package cromwell.backend.google.pipelines.v2beta.api

import cloud.nio.impl.drs.DrsConfig
import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import com.typesafe.config.ConfigFactory
import cromwell.backend.google.pipelines.common.action.ActionCommands.localizeFile
import cromwell.backend.google.pipelines.common.action.ActionLabels._
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.PipelinesApiFileInput
import cromwell.backend.google.pipelines.common.PipelinesApiJobPaths._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2beta.PipelinesConversions._
import cromwell.backend.google.pipelines.v2beta.ToParameter.ops._
import cromwell.backend.google.pipelines.v2beta.api.ActionBuilder.{EnhancedAction, cloudSdkShellAction}
import cromwell.core.path.Path
import cromwell.filesystems.drs.DrsPath

import scala.jdk.CollectionConverters._


trait Localization {
  def localizeActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount])
                     (implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Action] = {
    val localizationLabel = Map(Key.Tag -> Value.Localization)

    val gcsTransferLibraryContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsTransferLibraryName)
    val localizeGcsTransferLibrary = cloudSdkShellAction(localizeFile(
      cloudPath = createPipelineParameters.cloudCallRoot / GcsTransferLibraryName,
      containerPath = gcsTransferLibraryContainerPath))(mounts = mounts, labels = localizationLabel)

    val gcsLocalizationContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsLocalizationScriptName)
    val localizeGcsLocalizationScript = cloudSdkShellAction(localizeFile(
      cloudPath = createPipelineParameters.cloudCallRoot / GcsLocalizationScriptName,
      containerPath = gcsLocalizationContainerPath))(mounts = mounts, labels = localizationLabel)

    val gcsDelocalizationContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsDelocalizationScriptName)
    val localizeGcsDelocalizationScript = cloudSdkShellAction(localizeFile(
      cloudPath = createPipelineParameters.cloudCallRoot / GcsDelocalizationScriptName,
      containerPath = gcsDelocalizationContainerPath))(mounts = mounts, labels = localizationLabel)

    val runGcsLocalizationScript = cloudSdkShellAction(
      s"/bin/bash $gcsLocalizationContainerPath")(mounts = mounts, labels = localizationLabel)

    val drsInputs: List[DrsPath] = createPipelineParameters.inputOutputParameters.fileInputParameters.collect {
      case PipelinesApiFileInput(_, drsPath: DrsPath, _, _) => drsPath
    }

    val drsLocalizationActions = if (drsInputs.nonEmpty) {
      val drsLocalizationManifestContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(DrsLocalizationManifestName)
      val localizeDrsLocalizationManifest = cloudSdkShellAction(localizeFile(
        cloudPath = createPipelineParameters.cloudCallRoot / DrsLocalizationManifestName,
        containerPath = drsLocalizationManifestContainerPath))(mounts = mounts, labels = localizationLabel)

      // Requester pays project id is stored on each DrsPath, but will be the same for all DRS inputs to a
      // particular workflow because it's determined by the Google project set in workflow options.
      val requesterPaysProjectId: Option[String] = drsInputs.flatMap(_.requesterPaysProjectIdOption).headOption
      val runDrsLocalization = Localization.drsAction(drsLocalizationManifestContainerPath, mounts, localizationLabel, requesterPaysProjectId)
      List(localizeDrsLocalizationManifest, runDrsLocalization)
    } else List[Action]()

    // Any "classic" PAPI v2 one-at-a-time localizations for non-GCS inputs.
    val singletonLocalizations = createPipelineParameters.inputOutputParameters.fileInputParameters.flatMap(_.toActions(mounts).toList)

    val localizations =
      localizeGcsTransferLibrary ::
        localizeGcsLocalizationScript :: runGcsLocalizationScript ::
        drsLocalizationActions :::
        localizeGcsDelocalizationScript ::
        singletonLocalizations

    ActionBuilder.annotateTimestampedActions("localization", Value.Localization)(localizations)
  }
}

object Localization {

  def drsAction(manifestPath: Path,
                mounts: List[Mount],
                labels: Map[String, String],
                requesterPaysProjectId: Option[String]
               ): Action = {
    val config = ConfigFactory.load
    val drsResolverConfig = config.getConfig("filesystems.drs.global.config.resolver")
    val drsConfig = DrsConfig.fromConfig(drsResolverConfig)
    val drsDockerImage = config.getString("drs.localization.docker-image")

    val manifestArg = List("-m", manifestPath.pathAsString)
    val requesterPaysArg = requesterPaysProjectId.map(r => List("-r", r)).getOrElse(List.empty)
    val drsCommand = manifestArg ++ requesterPaysArg

    val drsResolverEnv = DrsConfig.toEnv(drsConfig)
    ActionBuilder
      .withImage(drsDockerImage)
      .withCommand(drsCommand: _*)
      .withMounts(mounts)
      .setEnvironment(drsResolverEnv.asJava)
      .withLabels(labels)
  }
}
