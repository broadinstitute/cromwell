package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.Mount
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common.PipelinesApiJobPaths._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.ToParameter.ops._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels.{Key, Value}
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.cloudSdkShellAction
import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands.localizeFile


trait Localization {
  def localizeActions(createPipelineParameters: CreatePipelineParameters, mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration) = {
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

    // Any "classic" PAPI v2 one-at-a-time localizations for non-GCS inputs.
    val singletonLocalizations = createPipelineParameters.inputOutputParameters.fileInputParameters.flatMap(_.toActions(mounts).toList)

    val localizations =
      localizeGcsTransferLibrary ::
        localizeGcsLocalizationScript :: runGcsLocalizationScript ::
        localizeGcsDelocalizationScript ::
        singletonLocalizations

    ActionBuilder.annotateTimestampedActions("localization", Value.Localization)(localizations)
  }
}
