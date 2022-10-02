package cromwell.backend.google.pipelines.common.action

import com.typesafe.config.ConfigFactory
import org.apache.commons.text.StringEscapeUtils
import net.ceedubs.ficus.Ficus._

object ActionUtils {
  /** Image to use for ssh access. */
  val sshImage = "gcr.io/cloud-genomics-pipelines/tools"

  /** Entry point on the ssh image. */
  val sshEntryPoint = "ssh-server"

  /** Port mappings for the ssh container. */
  val sshPortMappings: Map[String, Integer] = Map("22" -> Int.box(22))

  /***
   * A cloud-sdk 404.0.0-alpine + jq image has a decompressed size of 825 MB.
   *
   * {{{
   * % cat Dockerfile
   * FROM gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine
   * RUN apk update && apk upgrade && apk add --no-cache jq
   * % docker build . --quiet --tag cloud_sdk_plus_jq
   * sha256:8e6cc12a28ad1802253a593e6fdc5ec5b4cbf7d9b6445ca84f04547a81b40a92
   * % docker images | head -2
   * REPOSITORY                                  TAG                            IMAGE ID       CREATED          SIZE
   * cloud_sdk_plus_jq                           latest                         8e6cc12a28ad   8 seconds ago    825MB
   * }}}
   */
  val cromwellImagesSizeRoundedUpInGB = 1

  private val config = ConfigFactory.load().getConfig("google")

  /**
    * An image with the Google Cloud SDK installed.
    * http://gcr.io/google.com/cloudsdktool/cloud-sdk
    *
    * FYI additional older versions are available on DockerHub at:
    * https://hub.docker.com/r/google/cloud-sdk
    */
  val CloudSdkImage: String =
    config.getOrElse("cloud-sdk-image-url", "gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine")

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

  def timestampedMessage(message: String): String =
    s"""printf '%s %s\\n' "$$(date -u '+%Y/%m/%d %H:%M:%S')" ${shellEscaped(message)}"""

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
