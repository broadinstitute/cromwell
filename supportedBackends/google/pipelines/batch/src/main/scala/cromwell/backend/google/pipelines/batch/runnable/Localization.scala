package cromwell.backend.google.pipelines.batch.runnable

//import cloud.nio.impl.drs.DrsConfig
import com.google.cloud.batch.v1.{Runnable, Volume}
import com.typesafe.config.ConfigFactory
import cromwell.backend.google.pipelines.batch.api.GcpBatchRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.batch.models.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.batch.models.{GcpBatchFileInput, GcpBatchJobPaths}
import cromwell.backend.google.pipelines.batch.util.GcpBatchParameterConversions._
import cromwell.backend.google.pipelines.batch.util.ToParameter.ops._
import cromwell.core.path.Path
import cromwell.filesystems.drs.DrsPath
//import scala.jdk.CollectionConverters._

trait Localization {
  import GcpBatchJobPaths._
  import RunnableBuilder._
  import RunnableCommands._
  import RunnableLabels._

  def localizeRunnables(createPipelineParameters: CreatePipelineParameters, volumes: List[Volume])
                      (implicit gcsTransferConfiguration: GcsTransferConfiguration): List[Runnable] = {
    val localizationLabel = Map(Key.Tag -> Value.Localization)

    val gcsTransferLibraryContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsTransferLibraryName)
    val localizeGcsTransferLibrary = cloudSdkShellRunnable(localizeFile(
      cloudPath = createPipelineParameters.cloudCallRoot / GcsTransferLibraryName,
      containerPath = gcsTransferLibraryContainerPath))(volumes = volumes, labels = localizationLabel, flags = List.empty)

    val gcsLocalizationContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsLocalizationScriptName)
    val localizeGcsLocalizationScript = cloudSdkShellRunnable(localizeFile(
      cloudPath = createPipelineParameters.cloudCallRoot / GcsLocalizationScriptName,
      containerPath = gcsLocalizationContainerPath))(volumes = volumes, labels = localizationLabel, flags = List.empty)

    val gcsDelocalizationContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsDelocalizationScriptName)
    val localizeGcsDelocalizationScript = cloudSdkShellRunnable(localizeFile(
      cloudPath = createPipelineParameters.cloudCallRoot / GcsDelocalizationScriptName,
      containerPath = gcsDelocalizationContainerPath))(volumes = volumes, labels = localizationLabel, flags = List.empty)

    val runGcsLocalizationScript = cloudSdkShellRunnable(
      s"/bin/bash $gcsLocalizationContainerPath")(volumes = volumes, labels = localizationLabel, flags = List.empty)

    val drsInputs: List[DrsPath] = createPipelineParameters.inputOutputParameters.fileInputParameters.collect {
      case GcpBatchFileInput(_, drsPath: DrsPath, _, _) => drsPath
    }

    val drsLocalizationRunnable = if (drsInputs.nonEmpty) {
      val drsLocalizationManifestContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(DrsLocalizationManifestName)
      val localizeDrsLocalizationManifest = cloudSdkShellRunnable(localizeFile(
        cloudPath = createPipelineParameters.cloudCallRoot / DrsLocalizationManifestName,
        containerPath = drsLocalizationManifestContainerPath))(volumes = volumes, labels = localizationLabel, flags = List.empty)

      // Requester pays project id is stored on each DrsPath, but will be the same for all DRS inputs to a
      // particular workflow because it's determined by the Google project set in workflow options.
      val requesterPaysProjectId: Option[String] = drsInputs.flatMap(_.requesterPaysProjectIdOption).headOption
      val runDrsLocalization = Localization.drsRunnable(drsLocalizationManifestContainerPath, localizationLabel, requesterPaysProjectId)
      List(localizeDrsLocalizationManifest, runDrsLocalization)
    } else List[Runnable.Builder]()

    // Any "classic" PAPI v2 one-at-a-time localizations for non-GCS inputs.
    val singletonLocalizations = createPipelineParameters.inputOutputParameters.fileInputParameters.flatMap(_.toRunnables(volumes))

    val localizations =
      localizeGcsTransferLibrary ::
        localizeGcsLocalizationScript :: runGcsLocalizationScript ::
        drsLocalizationRunnable :::
        localizeGcsDelocalizationScript ::
        singletonLocalizations

    RunnableBuilder
      .annotateTimestampedRunnable("localization", Value.Localization)(volumes, localizations)
      .map(_.build)
  }
}

object Localization {

  // TODO: Avoid loading the global config because Cromwell already loaded it
  private lazy val config = ConfigFactory.load

  def drsRunnable(manifestPath: Path,
                labels: Map[String, String],
                requesterPaysProjectId: Option[String]
               ): Runnable.Builder = {
    import RunnableBuilder.EnhancedRunnableBuilder

//    val marthaConfig = config.getConfig("filesystems.drs.global.config.martha")
//    val drsConfig = DrsConfig.fromConfig(marthaConfig)
    val drsDockerImage = config.getString("drs.localization.docker-image")

    val manifestArg = List("-m", manifestPath.pathAsString)
    val requesterPaysArg = requesterPaysProjectId.map(r => List("-r", r)).getOrElse(List.empty)
    val drsCommand = manifestArg ++ requesterPaysArg

//    val marthaEnv = DrsConfig.toEnv(drsConfig)
    RunnableBuilder
      .withImage(drsDockerImage)
      .withCommand(drsCommand: _*)
//      .setEnvironment(marthaEnv)
//      .withLabels(labels)
  }
}
