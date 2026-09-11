package cromwell.backend.google.batch.runnable

import com.google.cloud.batch.v1.Runnable
import com.typesafe.config.ConfigFactory
import net.ceedubs.ficus.Ficus._
import org.apache.commons.text.StringEscapeUtils

object RunnableUtils {

  private val config = ConfigFactory.load().getConfig("google")

  val MountPoint: String = "/mnt/disks/cromwell_root"

  /**
    * An image with the Google Cloud SDK installed.
    * http://gcr.io/google.com/cloudsdktool/cloud-sdk
    *
    * Google deletes tags from this repository one year after they are published, so this pin needs periodic bumping.
    * When it expires, every runnable below fails to pull and jobs die before producing any output.
    */
  val CloudSdkImage: String =
    config.getOrElse("cloud-sdk-image-url", "gcr.io/google.com/cloudsdktool/cloud-sdk:583.0.0-alpine")

  /**
    * Batch sets `CLOUDSDK_PYTHON=/usr/bin/python3` on every runnable, but the alpine Cloud SDK image ships Python at
    * /usr/local/bin/python3 only, so `gsutil` there dies with exit 127 before writing anything to stderr. Blanking the
    * variable lets `gcloud` and `gsutil` locate their own interpreter, which works regardless of which image (or which
    * operator-supplied `cloud-sdk-image-url`) we end up running. (AN-601, XKCD-1987)
    */
  val CloudSdkEnvironment: Map[String, String] = Map("CLOUDSDK_PYTHON" -> "")

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
                    isBackground: Runnable => Boolean
  ): List[Runnable] = {
    val toBeSortedRunnables = localization ++ userRunnable ++ memoryRetryRunnable ++ deLocalization
    val sortedRunnables = toBeSortedRunnables.sortWith { case (runnable, _) =>
      isBackground(runnable)
    }

    containerSetup ++ monitoringSetup ++ checkpointingStart ++ sortedRunnables ++ checkpointingShutdown ++ monitoringShutdown
  }
}
