package cromwell.backend.google.pipelines.v2alpha1.api

import java.util.UUID

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import common.util.StringUtil._
import cromwell.backend.google.pipelines.common.PipelinesApiAsyncBackendJobExecutionActor
import cromwell.backend.google.pipelines.common.PipelinesApiConfigurationAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common.PipelinesApiJobPaths.GcsDelocalizationScriptName
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.RuntimeOutputMapping
import cromwell.backend.google.pipelines.v2alpha1.ToParameter.ops._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels.{Key, Value}
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands._
import cromwell.backend.google.pipelines.v2alpha1.api.Delocalization._
import cromwell.core.path.{DefaultPathBuilder, Path}
import wom.runtime.WomOutputRuntimeExtractor

import scala.collection.JavaConverters._
import scala.concurrent.duration._

object Delocalization {
  private val logsRoot = "/google/logs"
}

trait Delocalization {

  private def aggregatedLog = s"$logsRoot/output"

  private def delocalizeLogsAction(gcsLogPath: Path)(implicit localizationConfiguration: LocalizationConfiguration) = {
    cloudSdkShellAction(
    delocalizeDirectory(DefaultPathBuilder.build(logsRoot).get, gcsLogPath, PipelinesApiAsyncBackendJobExecutionActor.plainTextContentType)
    )(flags = List(ActionFlag.AlwaysRun), labels = Map(Key.Tag -> Value.Delocalization))
  }

  // Used for the final copy of the logs to make sure we have the most up to date version before terminating the job
  private def copyAggregatedLogToLegacyPath(callExecutionContainerRoot: Path, gcsLegacyLogPath: Path)(implicit localizationConfiguration: LocalizationConfiguration): Action = {
    cloudSdkShellAction(
      delocalizeFileTo(DefaultPathBuilder.build(aggregatedLog).get, gcsLegacyLogPath, PipelinesApiAsyncBackendJobExecutionActor.plainTextContentType)
    )(flags = List(ActionFlag.AlwaysRun), labels = Map(Key.Tag -> Value.Delocalization))
  }

  // Periodically copies the logs out to GCS
  private def copyAggregatedLogToLegacyPathPeriodic(callExecutionContainerRoot: Path, gcsLegacyLogPath: Path)(implicit localizationConfiguration: LocalizationConfiguration): Action = {
    cloudSdkShellAction(
      every(30.seconds) { delocalizeFileTo(DefaultPathBuilder.build(aggregatedLog).get, gcsLegacyLogPath, PipelinesApiAsyncBackendJobExecutionActor.plainTextContentType) }
    )(flags = List(ActionFlag.RunInBackground), labels = Map(Key.Tag -> Value.Background))
  }

  private def runtimeOutputExtractorAction(containerCallRoot: String,
                                    outputFile: String,
                                    mounts: List[Mount],
                                    womOutputRuntimeExtractor: WomOutputRuntimeExtractor): Action = {
    val commands = List(
      "-c",
      // Create the directory where the fofn will be written
      s"mkdir -p $$(dirname $outputFile) && " +
      s"cd $containerCallRoot && " +
      """echo "Runtime output files to be delocalized:" && """ +
      s"${womOutputRuntimeExtractor.command} | tee $outputFile"
    )

    ActionBuilder
      .withImage(womOutputRuntimeExtractor.dockerImage.getOrElse(cloudSdkImage))
      .setCommands(commands.asJava)
      .withMounts(mounts)
      .setEntrypoint("/bin/bash")
      // Saves us some time if something else fails before we get to run this action
      .withFlags(List(ActionFlag.DisableImagePrefetch))
      .withLabels(Map(Key.Tag -> Value.Delocalization))
  }

  private def delocalizeRuntimeOutputsScript(fofnPath: String, workflowRoot: Path, cloudCallRoot: Path)(implicit localizationConfiguration: LocalizationConfiguration) = {
    val gsutilCommand: String => String = { flag =>
      s"""rm -f $$HOME/.config/gcloud/gce && gsutil -m $flag cp -r $$line "${cloudCallRoot.pathAsString.ensureSlashed}$$gcs_path""""
    }

    def sedStripPrefix(prefix: String) = s"""sed -e "s/^${prefix.ensureSedEscaped}//""""

    // See RuntimeOutputMapping.prefixFilters for more info on why this is needed
    val prefixFilters = RuntimeOutputMapping
      .prefixFilters(workflowRoot)
      .map(sedStripPrefix)
      .mkString(" | ")

    /*
     * Delocalize all the files returned by the runtime output extractor
     */
    s"""|#!/bin/bash
        |
        |set -x
        |
        |if [ -f $fofnPath ]; then
        |  while IFS= read line
        |  do
        |    gcs_path=$$(echo $$line | $prefixFilters)
        |    (
        |       ${retry(recoverRequesterPaysError(cloudCallRoot)(gsutilCommand))}
        |    )
        |  done  <$fofnPath
        |fi""".stripMargin
  }

  private def delocalizeRuntimeOutputsAction(cloudCallRoot: Path, inputFile: String, workflowRoot: Path, mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration): Action = {
    val command = multiLineCommand(delocalizeRuntimeOutputsScript(inputFile, workflowRoot, cloudCallRoot))
    ActionBuilder.cloudSdkShellAction(command)(mounts, flags = List(ActionFlag.DisableImagePrefetch), labels = Map(Key.Tag -> Value.Delocalization))
  }

  def deLocalizeActions(createPipelineParameters: CreatePipelineParameters,
                        mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration): List[Action] = {
    val cloudCallRoot = createPipelineParameters.cloudCallRoot
    val callExecutionContainerRoot = createPipelineParameters.commandScriptContainerPath.parent

    val gcsLogDirectoryPath = createPipelineParameters.cloudCallRoot / "pipelines-logs"
    val gcsLegacyLogPath = createPipelineParameters.logGcsPath

    /*
     * Ideally temporaryFofnForRuntimeOutputFiles should be somewhere else than the execution directory (we could mount anther directory)
     * However because it runs after everything else there's no risk of polluting the task's results and the random ID ensures we don't override anything
     */
    val temporaryFofnDirectoryForRuntimeOutputFiles = callExecutionContainerRoot.pathAsString.ensureSlashed + UUID.randomUUID().toString.split("-")(0)
    val temporaryFofnForRuntimeOutputFiles = temporaryFofnDirectoryForRuntimeOutputFiles + "/runtime_output_files.txt"

    val runtimeExtractionActions = createPipelineParameters.womOutputRuntimeExtractor.toList flatMap { extractor =>
      List (
        runtimeOutputExtractorAction(callExecutionContainerRoot.pathAsString, temporaryFofnForRuntimeOutputFiles, mounts, extractor),
        delocalizeRuntimeOutputsAction(cloudCallRoot, temporaryFofnForRuntimeOutputFiles, createPipelineParameters.cloudWorkflowRoot, mounts)
      )
    }

    val gcsDelocalizationContainerPath = createPipelineParameters.commandScriptContainerPath.sibling(GcsDelocalizationScriptName)

    val delocalizationLabel = Map(Key.Tag -> Value.Delocalization)
    val runGcsDelocalizationScript: Action = cloudSdkShellAction(
      s"/bin/bash $gcsDelocalizationContainerPath")(mounts = mounts, labels = delocalizationLabel)

    ActionBuilder.annotateTimestampedActions("delocalization", Value.Delocalization)(
      runGcsDelocalizationScript ::
        createPipelineParameters.outputParameters.flatMap(_.toActions(mounts).toList) ++
        runtimeExtractionActions
    ) :+
      copyAggregatedLogToLegacyPath(callExecutionContainerRoot, gcsLegacyLogPath) :+
      copyAggregatedLogToLegacyPathPeriodic(callExecutionContainerRoot, gcsLegacyLogPath) :+
      delocalizeLogsAction(gcsLogDirectoryPath)
  }
}
