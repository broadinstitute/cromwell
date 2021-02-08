package cromwell.backend.google.pipelines.common.action

import org.apache.commons.text.StringEscapeUtils

object ActionUtils {
  /** Image to use for ssh access. */
  val sshImage = "gcr.io/cloud-genomics-pipelines/tools"

  /** Entry point on the ssh image. */
  val sshEntryPoint = "ssh-server"

  /** Port mappings for the ssh container. */
  val sshPortMappings = Map("22" -> Int.box(22))

  /*
   * At the moment, cloud-sdk (924MB for 276.0.0-slim) and stedolan/jq (182MB) decompressed ~= 1.1 GB
   */
  val cromwellImagesSizeRoundedUpInGB = 1

  /** Quotes a string such that it's compatible as a string argument in the shell. */
  def shellEscaped(any: Any): String = {
    val str = String.valueOf(any)
    /*
    NOTE: escapeXSI is more compact than wrapping in single quotes. Newlines are also stripped by the shell, as they
    are by escapeXSI. If for some reason escapeXSI doesn't 100% work, say because it ends up stripping some required
    newlines, then consider adding a check for newlines and then using:

    "'" + str.replace("'", "'\"'\"'") + "'"
     */
    StringEscapeUtils.escapeXSI(str)
  }

  /** Cloud Life Sciences / Genomics Pipelines API logs directory. */
  val logsRoot = "/google/logs"

  /** Cloud Life Sciences / Genomics Pipelines API logs file. */
  val aggregatedLog = s"$logsRoot/output"

  /**
    * Define a shared PID namespace for background action containers and their termination controller.
    * The value is "monitoring" for historical (first usage) reasons.
    */
  val backgroundActionPidNamespace = "monitoring"

  /**
    * monitoringTerminationAction is needed to gracefully terminate monitoring action,
    * because PAPIv2 currently sends SIGKILL to terminate background actions.
    *
    * A fixed timeout is used to avoid hunting for monitoring PID.
    */
  private val backgroundActionTerminationGraceTime = 10

  val terminateAllBackgroundActionsCommand: String = s"kill -TERM -1 && sleep $backgroundActionTerminationGraceTime || true"

  /** Start background actions first, leave the rest as is */
  def sortActions[Action](containerSetup: List[Action],
                          localization: List[Action],
                          userAction: List[Action],
                          memoryRetryAction: List[Action],
                          deLocalization: List[Action],
                          monitoringSetup: List[Action],
                          monitoringShutdown: List[Action],
                          checkpointingStart: List[Action],
                          checkpointingShutdown: List[Action],
                          sshAccess: List[Action],
                          isBackground: Action => Boolean,
                         ): List[Action] = {
    val toBeSortedActions = localization ++ userAction ++ memoryRetryAction ++ deLocalization
    val sortedActions = toBeSortedActions.sortWith({
      case (action, _) => isBackground(action)
    })

    sshAccess ++ containerSetup ++ monitoringSetup ++ checkpointingStart ++ sortedActions ++ checkpointingShutdown ++ monitoringShutdown
  }
}
