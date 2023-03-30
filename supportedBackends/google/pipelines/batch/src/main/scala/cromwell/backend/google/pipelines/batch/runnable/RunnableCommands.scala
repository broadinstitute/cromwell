package cromwell.backend.google.pipelines.batch.runnable

//import java.util.UUID
//import akka.http.scaladsl.model.ContentType
//import common.util.StringUtil._

import scala.concurrent.duration._
import mouse.all._
import cromwell.core.path.Path
import cromwell.filesystems.gcs.GcsPath
import org.apache.commons.text.StringEscapeUtils
import cromwell.backend.google.pipelines.batch.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.pipelines.common.action.ActionUtils._

import cromwell.filesystems.gcs.RequesterPaysErrors._

//import org.apache.commons.codec.binary.Base64
//import org.apache.commons.text.StringEscapeUtils

//import java.nio.charset.StandardCharsets
//import scala.concurrent.duration._

object RunnableCommands {
  implicit val waitBetweenRetries: FiniteDuration = 5.seconds

  implicit class EnhancedCromwellPath(val path: Path) extends AnyVal {
    def projectId: String = path match {
      case gcs: GcsPath => gcs.projectId
      case _ => ""
    }
  }

  implicit class ShellPath(val path: Path) extends AnyVal {
    // The command String runs in Bourne shell so shell metacharacters in filenames must be escaped
    def escape: String = StringEscapeUtils.escapeXSI(path.pathAsString)
  }

  def retry(f: => String)(implicit gcsTransferConfiguration: GcsTransferConfiguration, wait: FiniteDuration): String = {
    s"""for i in $$(seq ${gcsTransferConfiguration.transferAttempts}); do
       |  (
       |    $f
       |  )
       |  RC=$$?
       |  if [ "$$RC" = "0" ]; then
       |    break
       |  fi
       |  if [ $$i -lt ${gcsTransferConfiguration.transferAttempts} ]; then
       |    ${s"""Waiting ${wait.toSeconds} seconds and retrying""" |> timestampedMessage}
       |    sleep ${wait.toSeconds}
       |  fi
       |done
       |exit "$$RC"""".stripMargin
  }

  def localizeFile(cloudPath: Path, containerPath: Path)
                  (implicit gcsTransferConfiguration: GcsTransferConfiguration): String = {
    retry {
      recoverRequesterPaysError(cloudPath) { flag =>
        s"rm -f $$HOME/.config/gcloud/gce && gsutil $flag cp ${cloudPath.escape} ${containerPath.escape}"
      }
    }
  }

  def recoverRequesterPaysError(path: Path)(f: String => String): String = {
    val commandWithoutProject = f("")
    val commandWithProject = f(s"-u ${path.projectId}")

    s"""$commandWithoutProject > gsutil_output.txt 2>&1
       |# Record the exit code of the gsutil command without project flag
       |RC_GSUTIL=$$?
       |if [ "$$RC_GSUTIL" != "0" ]; then
       |  ${s"$commandWithoutProject failed" |> timestampedMessage}
       |  # Print the reason of the failure
       |  cat gsutil_output.txt
       |
       |  # Check if it matches the BucketIsRequesterPaysErrorMessage
       |  if grep -q "$BucketIsRequesterPaysErrorMessage" gsutil_output.txt; then
       |    ${"Retrying with user project" |> timestampedMessage}
       |    $commandWithProject
       |  else
       |    exit "$$RC_GSUTIL"
       |  fi
       |else
       |  exit 0
       |fi""".stripMargin
  }

}
