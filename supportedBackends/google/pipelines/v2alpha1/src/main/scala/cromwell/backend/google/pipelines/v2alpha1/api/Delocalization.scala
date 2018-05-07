package cromwell.backend.google.pipelines.v2alpha1.api

import com.google.api.services.genomics.v2alpha1.model.{Action, Mount}
import cromwell.backend.google.pipelines.common.api.PipelinesApiRequestFactory.CreatePipelineParameters
import cromwell.backend.google.pipelines.v2alpha1.PipelinesConversions._
import cromwell.backend.google.pipelines.v2alpha1.api.ActionBuilder._
import cromwell.backend.google.pipelines.v2alpha1.api.Delocalization._

object Delocalization {
  private val logsRoot = "/google/logs"
  private val actionsLogRoot = logsRoot + "/action"
}

trait Delocalization {
  private def actionLogRoot(number: Int) = s"$actionsLogRoot/$number"

  private def stdout(number: Int) = s"${actionLogRoot(number)}/stdout"
  private def stderr(number: Int) = s"${actionLogRoot(number)}/stderr"
  private def aggregatedLog = s"$logsRoot/output"

  private def delocalizeLogsAction(gcsLogPath: String) = {
    gsutilAsText("cp", "-r", "/google/logs", gcsLogPath)(flags = List(ActionFlag.AlwaysRun))
  }

  // The logs are now located in the pipelines-logs directory
  // To keep the behavior similar to V1, we copy stdout/stderr from the user action to the call directory,
  // along with the aggregated log file
  private def copyLogsToLegacyPaths(stdoutPath: String, stderrPath: String, userActionNumber: Int, gcsLegacyLogPath: String) = List (
    gsutilAsText("cp", stdout(userActionNumber), stdoutPath)(flags = List(ActionFlag.AlwaysRun)),
    gsutilAsText("cp", stderr(userActionNumber), stderrPath)(flags = List(ActionFlag.AlwaysRun)),
    gsutilAsText("cp", aggregatedLog, gcsLegacyLogPath)(flags = List(ActionFlag.AlwaysRun))
  )

  /**
    * The user action number is the index of the user's action in the list of actions
    * It's used to copy stdout / stderr and the logs back to the execution directory
    */
  def deLocalizeActions(createPipelineParameters: CreatePipelineParameters,
                        mounts: List[Mount],
                        userActionNumber: Int): List[Action] = {
    val gcsLogDirectoryPath = createPipelineParameters.cloudCallRoot / "pipelines-logs"
    val gcsLegacyLogPath = createPipelineParameters.logGcsPath.pathAsString

    val stdoutPath = createPipelineParameters.logGcsPath
      .sibling(createPipelineParameters.logGcsPath.nameWithoutExtensionNoIo + "-stdout.log")
      .pathAsString

    val stderrPath = createPipelineParameters.logGcsPath
      .sibling(createPipelineParameters.logGcsPath.nameWithoutExtensionNoIo + "-stderr.log")
      .pathAsString

    createPipelineParameters.outputParameters.map(_.toAction(mounts)) ++
      copyLogsToLegacyPaths(stdoutPath, stderrPath, userActionNumber, gcsLegacyLogPath) :+
      delocalizeLogsAction(gcsLogDirectoryPath.pathAsString)
  }
}
