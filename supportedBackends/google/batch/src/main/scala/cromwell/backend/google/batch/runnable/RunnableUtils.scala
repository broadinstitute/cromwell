package cromwell.backend.google.batch.runnable

import com.google.cloud.batch.v1.Runnable
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._
import org.apache.commons.text.StringEscapeUtils

object RunnableUtils {

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

  private val config = ConfigFactory.load().getConfig("google")

  /**
    * An image with the Google Cloud SDK installed.
    * http://gcr.io/google.com/cloudsdktool/cloud-sdk
    *
    * FYI additional older versions are available on DockerHub at:
    * https://hub.docker.com/r/google/cloud-sdk
    *
    * When updating this value, also consider updating the CromwellImagesSizeRoundedUpInGB below.
    */
  val CloudSdkImage: String =
    // config.getOrElse("cloud-sdk-image-url", "gcr.io/google.com/cloudsdktool/cloud-sdk:354.0.0-alpine")
    config.getOrElse("cloud-sdk-image-url", "gcr.io/google.com/cloudsdktool/cloud-sdk:434.0.0-alpine")
  /*
   * At the moment, cloud-sdk (584MB for 354.0.0-alpine) and stedolan/jq (182MB) decompressed ~= 0.8 GB
   */
  val CromwellImagesSizeRoundedUpInGB: Int =
    config.getOrElse("cloud-sdk-image-size-gb", 1)

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

  private val backgroundRunnableTerminationGraceTime = 10

  val terminateAllBackgroundRunnablesCommand: String =
    s"kill -TERM -1 && sleep $backgroundRunnableTerminationGraceTime || true"

  def timestampedMessage(message: String): String =
    s"""printf '%s %s\\n' "$$(date -u '+%Y/%m/%d %H:%M:%S')" ${shellEscaped(message)}"""

  /** Start background runnables first, leave the rest as is */
  def sortRunnables(containerSetup: List[Runnable],
                    localization: List[Runnable],
                    userRunnable: List[Runnable],
                    memoryRetryRunnable: List[Runnable],
                    deLocalization: List[Runnable],
                    monitoringSetup: List[Runnable],
                    monitoringShutdown: List[Runnable],
                    checkpointingStart: List[Runnable],
                    checkpointingShutdown: List[Runnable],
                    sshAccess: List[Runnable],
                    isBackground: Runnable => Boolean
  ): List[Runnable] = {
    val toBeSortedRunnables = localization ++ userRunnable ++ memoryRetryRunnable ++ deLocalization
    val sortedRunnables = toBeSortedRunnables.sortWith { case (runnable, _) =>
      isBackground(runnable)
    }

    sshAccess ++ containerSetup ++ monitoringSetup ++ checkpointingStart ++ sortedRunnables ++ checkpointingShutdown ++ monitoringShutdown
  }
}
