package cromwell.backend.google.pipelines.v2beta.api

import com.google.api.services.lifesciences.v2beta.model.{Action, Mount}
import com.typesafe.config.ConfigFactory
import cromwell.backend.google.pipelines.common.action.ActionCommands.localizeFile
import cromwell.backend.google.pipelines.common.action.ActionLabels._
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.PipelinesApiJobPaths._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2beta.PipelinesConversions._
import cromwell.backend.google.pipelines.v2beta.ToParameter.ops._
import cromwell.backend.google.pipelines.v2beta.api.ActionBuilder.{EnhancedAction, cloudSdkShellAction}
import cromwell.core.path.Path


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

    val drsLocalizationManifestContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(DrsLocalizationManifestName)
    val localizeDrsLocalizationManifest = cloudSdkShellAction(localizeFile(
      cloudPath = createPipelineParameters.cloudCallRoot / DrsLocalizationManifestName,
      containerPath = drsLocalizationManifestContainerPath))(mounts = mounts, labels = localizationLabel)

    val runGcsLocalizationScript = cloudSdkShellAction(
      s"/bin/bash $gcsLocalizationContainerPath")(mounts = mounts, labels = localizationLabel)

    val runDrsLocalization = drsAction(createPipelineParameters, drsLocalizationManifestContainerPath, mounts, localizationLabel)

    // Any "classic" PAPI v2 one-at-a-time localizations for non-GCS inputs.
    val singletonLocalizations = createPipelineParameters.inputOutputParameters.fileInputParameters.flatMap(_.toActions(mounts).toList)

    val localizations =
      localizeGcsTransferLibrary ::
        localizeGcsLocalizationScript :: runGcsLocalizationScript ::
        localizeDrsLocalizationManifest :: runDrsLocalization ::
        localizeGcsDelocalizationScript ::
        singletonLocalizations

    ActionBuilder.annotateTimestampedActions("localization", Value.Localization)(localizations)
  }

  private def drsAction(createPipelineParameters: CreatePipelineParameters, manifestPath: Path, mounts: List[Mount], labels: Map[String, String]) = {
    lazy val config = ConfigFactory.load
    import scala.jdk.CollectionConverters._

    val drsDockerImage = config.getString("drs.localization.docker-image")
    val drsCommand = List("-m", manifestPath.pathAsString)
    // TODO: find a way to read this config without having a DrsPath on hand
    val marthaEnv = Map(
      "MARTHA_URL" -> "https://us-central1-broad-dsde-dev.cloudfunctions.net/martha_v3",
      "MARTHA_NUM_RETRIES" -> "3",
      "MARTHA_WAIT_INITIAL_SECONDS" -> "10",
      "MARTHA_WAIT_MAXIMUM_SECONDS" -> "30",
      "MARTHA_WAIT_MULTIPLIER" -> "1.5",
      "MARTHA_WAIT_RANDOMIZATION_FACTOR" -> "0.1"
    )
    ActionBuilder
      .withImage(drsDockerImage)
      .withCommand(drsCommand: _*)
      .withMounts(mounts)
      .setEnvironment(marthaEnv.asJava)
      .withLabels(labels)
  }
}
