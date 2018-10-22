package cromwell.backend.google.pipelines.v2alpha1.api

import java.util.UUID

import akka.http.scaladsl.model.ContentTypes
import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import common.util.StringUtil._
import cromwell.backend.google.pipelines.common.PipelinesApiAttributes.LocalizationConfiguration
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.ToParameter.ops._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder.Labels.Value
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionCommands._
import cromwell.backend.google.pipelines.v2alpha1.api.Delocalization._
import cromwell.core.path.{DefaultPathBuilder, Path}

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

  private def parseOutputJsonAction(containerCallRoot: String, outputDirectory: String, outputFile: String, mounts: List[Mount]): Action = {
    val commands = List(
      "-c",
      // Create the directory where the output of jq will be written
      s"mkdir -p $outputDirectory &&" +
      // Read the content of cwl.output.json - redirect stderr to /dev/null so it doesn't fail if there's no cwl.output.json
      s" cat ${containerCallRoot.ensureSlashed}cwl.output.json 2>/dev/null" +
        // Pipe the result to jq. Parse the json and traverse it looking for "path" or "location" values,
        // stripping any fully qualified file:// URLs.
      """ | jq -r '.. | .path? // .location? // empty | gsub("file://"; "")' """ +
      // Redirect the result to a file that can be used in the next action to delocalize the extracted files / directories
      s" > $outputFile"
    )

    ActionBuilder
      .withImage("stedolan/jq@sha256:a61ed0bca213081b64be94c5e1b402ea58bc549f457c2682a86704dd55231e09")
      .setCommands(commands.asJava)
      .withMounts(mounts)
      .setEntrypoint("/bin/bash")
      // Saves us some time if something else fails before we get to run this action
      .withFlags(List(ActionFlag.DisableImagePrefetch))
  }

  private def delocalizeOutputJsonFilesAction(cloudCallRoot: Path, inputFile: String, workflowRoot: String, mounts: List[Mount]): Action = {
    val sedStripSlashPrefix = "s/^\\///"
    val gsutilCommand: String => String = flag => s"""gsutil -m $flag cp -r % ${cloudCallRoot.pathAsString.ensureSlashed}$$(echo % | sed -e "$sedStripSlashPrefix")"""
    val command =
        // Read the file containing files to delocalize
        s"cat $inputFile 2>/dev/null" +
        /*
         * Pipe the result to xargs and execute a gsutil cp command, appending the path to the cloudCallRoot
         * Use a subshell so that we can strip the potential leading slash from the local path and avoid double slashes in the cloud path
         * It also seems that the gsutil fails without sub shelling
         * We can't use the gsutil -I flag here because it would lose the directory structure once it gets copied to the bucket
         * sh -c 'gsutil cp % $(echo % | sed -e "s/^\///")'
         */
        s""" | xargs -I % sh -c '${recoverRequesterPaysError(cloudCallRoot)(gsutilCommand)}'"""

    ActionBuilder.cloudSdkShellAction(command)(mounts)
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
    val parseAction = parseOutputJsonAction(callExecutionContainerRoot.pathAsString, temporaryFofnDirectoryForCwlOutputJson, temporaryFofnForCwlOutputJson, mounts)
    val delocalizeAction = delocalizeOutputJsonFilesAction(cloudCallRoot, temporaryFofnForCwlOutputJson, createPipelineParameters.cloudWorkflowRoot.pathAsString, mounts)

    ActionBuilder.annotateTimestampedActions("delocalization", Value.Delocalization)(
      createPipelineParameters.outputParameters.flatMap(_.toActions(mounts).toList) ++
        List(parseAction, delocalizeAction)
    ) :+
      copyAggregatedLogToLegacyPath(callExecutionContainerRoot, gcsLegacyLogPath) :+
      copyAggregatedLogToLegacyPathPeriodic(callExecutionContainerRoot, createPipelineParameters.logGcsPath) :+
      delocalizeLogsAction(gcsLogDirectoryPath)
  }
}
