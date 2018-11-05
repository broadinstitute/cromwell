package cromwell.backend.google.pipelines.v2alpha1.api

import java.util.UUID

import akka.http.scaladsl.model.ContentTypes
import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import common.util.StringUtil._
import cromwell.backend.google.pipelines.common.PipelinesApiAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.RuntimeOutputMapping
import cromwell.backend.google.pipelines.v2alpha1.ToParameter.ops._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels.Value
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands._
import cromwell.backend.google.pipelines.v2alpha1.api.Delocalization._
import cromwell.core.path.{DefaultPathBuilder, Path}
import org.apache.commons.codec.binary.Base64
import wom.runtime.WomOutputRuntimeExtractor

import scala.collection.JavaConverters._
import scala.concurrent.duration._

object Delocalization {
  private val logsRoot = "/google/logs"
  val plainTextContentType = Option(ContentTypes.`text/plain(UTF-8)`)
}

trait Delocalization {

  private def aggregatedLog = s"$logsRoot/output"

  private def delocalizeLogsAction(gcsLogPath: Path)(implicit localizationConfiguration: LocalizationConfiguration) = {
    cloudSdkShellAction(
    delocalizeDirectory(DefaultPathBuilder.build(logsRoot).get, gcsLogPath, plainTextContentType)
    )(flags = List(ActionFlag.AlwaysRun))
  }

  // Used for the final copy of the logs to make sure we have the most up to date version before terminating the job
  private def copyAggregatedLogToLegacyPath(callExecutionContainerRoot: Path, gcsLegacyLogPath: Path)(implicit localizationConfiguration: LocalizationConfiguration): Action = {
    cloudSdkShellAction(
      delocalizeFileTo(DefaultPathBuilder.build(aggregatedLog).get, gcsLegacyLogPath, plainTextContentType)
    )(flags = List(ActionFlag.AlwaysRun))
  }

  // Periodically copies the logs out to GCS
  private def copyAggregatedLogToLegacyPathPeriodic(callExecutionContainerRoot: Path, gcsLegacyLogPath: Path)(implicit localizationConfiguration: LocalizationConfiguration): Action = {
    cloudSdkShellAction(
      every(30.seconds) { delocalizeFileTo(DefaultPathBuilder.build(aggregatedLog).get, gcsLegacyLogPath, plainTextContentType) }
    )(flags = List(ActionFlag.RunInBackground))
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
      s"${womOutputRuntimeExtractor.command} > $outputFile"
    )

    ActionBuilder
      .withImage(womOutputRuntimeExtractor.dockerImage.getOrElse(cloudSdkImage))
      .setCommands(commands.asJava)
      .withMounts(mounts)
      .setEntrypoint("/bin/bash")
      // Saves us some time if something else fails before we get to run this action
      .withFlags(List(ActionFlag.DisableImagePrefetch))
  }
  
  private def delocalizeRuntimeOutputsScript(fofnPath: String, workflowRoot: Path, cloudCallRoot: Path) = {
    val gsutilCommand: String => String = { flag =>
      s"""gsutil -m $flag cp -r $$line "${cloudCallRoot.pathAsString.ensureSlashed}$$gcs_path""""
    }
    
    def sedStripPrefix(prefix: String) = s"""sed -e "s/^${prefix.sedEscaped}//""""
    
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
        |    ${recoverRequesterPaysError(cloudCallRoot)(gsutilCommand)}
        |  done  <$fofnPath
        |fi""".stripMargin
  }

  private def delocalizeRuntimeOutputsAction(cloudCallRoot: Path, inputFile: String, workflowRoot: Path, mounts: List[Mount]): Action = {
    def multiLineCommand(commandString: String) = {
      val randomUuid = UUID.randomUUID().toString
      val base64EncodedScript = Base64.encodeBase64String(commandString.getBytes)
      val scriptPath = s"/tmp/$randomUuid.sh"

      s"apt-get install --assume-yes coreutils && " +
        s"echo $base64EncodedScript | base64 --decode > $scriptPath && " +
        s"chmod u+x $scriptPath && " +
        s"sh $scriptPath"
    }
    
    val command = multiLineCommand(delocalizeRuntimeOutputsScript(inputFile, workflowRoot, cloudCallRoot))
    ActionBuilder.cloudSdkShellAction(command)(mounts, flags = List(ActionFlag.DisableImagePrefetch))
  }

  def deLocalizeActions(createPipelineParameters: CreatePipelineParameters,
                        mounts: List[Mount])(implicit localizationConfiguration: LocalizationConfiguration): List[Action] = {
    val cloudCallRoot = createPipelineParameters.cloudCallRoot
    val callExecutionContainerRoot = createPipelineParameters.commandScriptContainerPath.parent

    val gcsLogDirectoryPath = createPipelineParameters.cloudCallRoot / "pipelines-logs"
    val gcsLegacyLogPath = createPipelineParameters.logGcsPath

    /*
     * CWL specific delocalization. For now this always runs, even for WDL jobs.
     * Ideally temporaryFofnForCwlOutputJson should be somewhere else than the execution directory (we could mount anther directory)
     * However because it runs after everything else there's no risk of polluting the task's results and the random ID ensures we don't override anything
     */
    val temporaryFofnDirectoryForCwlOutputJson = callExecutionContainerRoot.pathAsString.ensureSlashed + UUID.randomUUID().toString.split("-")(0)
    val temporaryFofnForCwlOutputJson = temporaryFofnDirectoryForCwlOutputJson + "/cwl_output_json_references.txt"

    val runtimeExtractionActions = createPipelineParameters.womOutputRuntimeExtractor.toList flatMap { extractor =>
      List (
        runtimeOutputExtractorAction(callExecutionContainerRoot.pathAsString, temporaryFofnForCwlOutputJson, mounts, extractor),
        delocalizeRuntimeOutputsAction(cloudCallRoot, temporaryFofnForCwlOutputJson, createPipelineParameters.cloudWorkflowRoot, mounts)
      )
    }

    ActionBuilder.annotateTimestampedActions("delocalization", Value.Delocalization)(
      createPipelineParameters.outputParameters.flatMap(_.toActions(mounts).toList) ++
        runtimeExtractionActions
    ) :+
      copyAggregatedLogToLegacyPath(callExecutionContainerRoot, gcsLegacyLogPath) :+
      copyAggregatedLogToLegacyPathPeriodic(callExecutionContainerRoot, createPipelineParameters.logGcsPath) :+
      delocalizeLogsAction(gcsLogDirectoryPath)
  }
}
