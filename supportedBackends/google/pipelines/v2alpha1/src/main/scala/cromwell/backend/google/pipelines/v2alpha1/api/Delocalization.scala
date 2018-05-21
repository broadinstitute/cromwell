package cromwell.backend.google.pipelines.v2alpha1.api

import java.util.UUID

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import common.util.StringUtil._
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.Delocalization._
import cromwell.backend.io.JobPaths
import cromwell.core.StandardPaths
import cromwell.core.path.Path

import scala.collection.JavaConverters._

object Delocalization {
  private val logsRoot = "/google/logs"
  private val actionsLogRoot = logsRoot + "/action"
}

trait Delocalization {
  private def actionLogRoot(number: Int) = s"$actionsLogRoot/$number"

  private def stdout(number: Int) = s"${actionLogRoot(number)}/stdout"
  private def stderr(number: Int) = s"${actionLogRoot(number)}/stderr"
  private def aggregatedLog = s"$logsRoot/output"

  private def delocalizeLogsAction(gcsLogPath: String, projectId: String) = {
    gsutilAsText("-u", projectId, "-m", "cp", "-r", "/google/logs", gcsLogPath)(flags = List(ActionFlag.AlwaysRun))
  }

  // The logs are now located in the pipelines-logs directory
  // To keep the behavior similar to V1, we copy stdout/stderr from the user action to the call directory,
  // along with the aggregated log file
  private def copyLogsToLegacyPaths(stdoutPath: String, stderrPath: String, userActionNumber: Int, gcsLegacyLogPath: String, projectId: String) = List (
    gsutilAsText("-u", projectId, "cp", stdout(userActionNumber), stdoutPath)(flags = List(ActionFlag.AlwaysRun)),
    gsutilAsText("-u", projectId, "cp", stderr(userActionNumber), stderrPath)(flags = List(ActionFlag.AlwaysRun)),
    gsutilAsText("-u", projectId, "cp", aggregatedLog, gcsLegacyLogPath)(flags = List(ActionFlag.AlwaysRun))
  )

  private def parseOutputJsonAction(containerCallRoot: String, outputDirectory: String, outputFile: String, mounts: List[Mount]): Action = {
    val commands = List(
      "-c",
      // Create the directory where the output of jq will be written
      s"mkdir -p $outputDirectory &&" +
      // Read the content of cwl.output.json - redirect stderr to /dev/null so it doesn't fail if there's no cwl.output.json
      s" cat ${containerCallRoot.ensureSlashed}cwl.output.json 2>/dev/null" +
        // Pipe the result to jq. Parse the json and traverse it looking for "path" values
      " | jq -r '.. | .path? // empty'" +
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

  private def delocalizeOutputJsonFilesAction(cloudCallRoot: String, inputFile: String, mounts: List[Mount]): Action = {
    val sedStripSlashPrefix = "s/^\\///"
    val commands = List(
      "/bin/bash",
      "-c",
        // Read the file containing files to delocalize
        s"cat $inputFile 2>/dev/null" +
        /*
         * Pipe the result to xargs and execute a gsutil cp command, appending the path to the cloudCallRoot
         * Use a subshell so that we can strip the potential leading slash from the local path and avoid double slashes in the cloud path
         * It also seems that the gsutil fails without sub shelling
         * We can't use the gsutil -I flag here because it would lose the directory structure once it gets copied to the bucket
         * sh -c 'gsutil cp % $(echo % | sed -e "s/^\///")'
         */
        s""" | xargs -I % sh -c 'gsutil -m cp -r % ${cloudCallRoot.ensureSlashed}$$(echo % | sed -e "$sedStripSlashPrefix")'"""
    )

    ActionBuilder
      .cloudSdkAction
      .setCommands(commands.asJava)
      .withMounts(mounts)
  }

  /**
    * The user action number is the index of the user's action in the list of actions
    * It's used to copy stdout / stderr and the logs back to the execution directory
    */
  def deLocalizeActions(createPipelineParameters: CreatePipelineParameters,
                        mounts: List[Mount],
                        userActionNumber: Int,
                        jobPaths: Option[JobPaths]): List[Action] = {
    val cloudCallRoot = createPipelineParameters.cloudCallRoot.pathAsString
    val callExecutionContainerRoot = createPipelineParameters.commandScriptContainerPath.parent.pathAsString

    val gcsLogDirectoryPath = createPipelineParameters.cloudCallRoot / "pipelines-logs"
    val gcsLegacyLogPath = createPipelineParameters.logGcsPath.pathAsString

    val parent = createPipelineParameters.logGcsPath.parent

    val standardFile: (StandardPaths => Path, String) => String = (get, name) =>
      jobPaths.map(j => parent.resolve(get(j.standardPaths).getFileName))
        .getOrElse(parent.resolve(createPipelineParameters.logGcsPath.nameWithoutExtensionNoIo + s"-$name.log"))
        .pathAsString

    val List(stdoutPath, stderrPath) = List[(StandardPaths => Path, String)]((_.output, "stdout"), (_.error, "stderr")) map standardFile.tupled

    /*
     * CWL specific delocalization. For now this always runs, even for WDL jobs.
     * Ideally temporaryFofnForCwlOutputJson should be somewhere else than the execution directory (we could mount anther directory)
     * However because it runs after everything else there's no risk of polluting the task's results and the random ID ensures we don't override anything
     */
    val temporaryFofnDirectoryForCwlOutputJson = callExecutionContainerRoot.ensureSlashed + UUID.randomUUID().toString.split("-")(0)
    val temporaryFofnForCwlOutputJson = temporaryFofnDirectoryForCwlOutputJson + "/cwl_output_json_references.txt"
    val parseAction = parseOutputJsonAction(callExecutionContainerRoot, temporaryFofnDirectoryForCwlOutputJson, temporaryFofnForCwlOutputJson, mounts)
    val delocalizeAction = delocalizeOutputJsonFilesAction(cloudCallRoot, temporaryFofnForCwlOutputJson, mounts)

    val projectId = createPipelineParameters.projectId

    createPipelineParameters.outputParameters.map(_.toAction(mounts, projectId)) ++
      List(parseAction, delocalizeAction) ++
      copyLogsToLegacyPaths(stdoutPath, stderrPath, userActionNumber, gcsLegacyLogPath, projectId) :+
      delocalizeLogsAction(gcsLogDirectoryPath.pathAsString, projectId)
  }
}
