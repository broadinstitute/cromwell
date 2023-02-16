package drs.localizer.downloaders
import cats.effect.{ExitCode, IO}
import com.typesafe.scalalogging.StrictLogging
import drs.localizer.downloaders.GcsUriDownloader.RequesterPaysErrorMsg

import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Path}
import scala.sys.process.{Process, ProcessLogger}

case class GcsUriDownloader(gcsUrl: String,
                            serviceAccountJson: Option[String],
                            downloadLoc: String,
                            requesterPaysProjectIdOption: Option[String]) extends Downloader with StrictLogging {

  override def download: IO[DownloadResult] = {

    logger.info(s"Requester Pays project ID is $requesterPaysProjectIdOption")
    logger.info(s"Attempting to download $gcsUrl to $downloadLoc")

    val copyProcess = serviceAccountJson match {
      case Some(sa) =>
        // if DRS Resolver returned a SA, use that SA for gsutil instead of default credentials
        val tempCredentialDir: Path = Files.createTempDirectory("gcloudTemp_").toAbsolutePath
        val saJsonPath: Path = tempCredentialDir.resolve("sa.json")
        Files.write(saJsonPath, sa.getBytes(StandardCharsets.UTF_8))
        val extraEnv = Map("CLOUDSDK_CONFIG" -> tempCredentialDir.toString)
        val copyCommand = Seq("bash", "-c", generateDownloadScript(gcsUrl, Option(saJsonPath)))
        Process(copyCommand, None, extraEnv.toSeq: _*)
      case None =>
        // No SA returned from DRS Resolver. gsutil will use the application default credentials.
        val copyCommand = Seq("bash", "-c", generateDownloadScript(gcsUrl, None))
        Process(copyCommand)
    }

    // run the multiple bash script to download file and log stream sent to stdout and stderr using ProcessLogger
    val returnCode = copyProcess ! ProcessLogger(logger.underlying.info, logger.underlying.error)

    val result = if (returnCode == 0) DownloadSuccess else RecognizedRetryableDownloadFailure(exitCode = ExitCode(returnCode))

    IO.pure(result)
  }

  /**
    * Bash to download the GCS file using `gsutil`.
    */
  def generateDownloadScript(gcsUrl: String, saJsonPathOption: Option[Path]): String = {

    def gcsCopyCommand(flag: String = ""): String = s"gsutil $flag cp $gcsUrl $downloadLoc"

    def setServiceAccount(): String = {
      saJsonPathOption match {
        case Some(saJsonPath) =>
          s"""# Set gsutil to use the service account returned from the DRS Resolver
             |gcloud auth activate-service-account --key-file=$saJsonPath > gcloud_output.txt 2>&1
             |RC_GCLOUD=$$?
             |if [ "$$RC_GCLOUD" != "0" ]; then
             |  echo "Failed to activate service account returned from the DRS Resolver. File won't be downloaded. Error: $$(cat gcloud_output.txt)" >&2
             |  exit "$$RC_GCLOUD"
             |else
             |  echo "Successfully activated service account; Will continue with download. $$(cat gcloud_output.txt)"
             |fi
             |""".stripMargin
        case None => ""
      }
    }

    def recoverWithRequesterPays(): String = {
      requesterPaysProjectIdOption match {
        case Some(userProject) =>
          s"""if [ "$$RC_GSUTIL" != "0" ]; then
             |  # Check if error is requester pays. If yes, retry gsutil copy using project flag
             |  if grep -q '$RequesterPaysErrorMsg' gsutil_output.txt; then
             |    echo "Received 'Bucket is requester pays' error. Attempting again using Requester Pays billing project"
             |    ${gcsCopyCommand(s"-u $userProject")} > gsutil_output.txt 2>&1
             |    RC_GSUTIL=$$?
             |  fi
             |fi
             |""".stripMargin
        case None => ""
      }
    }

    // bash to download the GCS file using gsutil
    s"""set -euo pipefail
       |set +e
       |
       |${setServiceAccount()}
       |
       |# Run gsutil copy without using project flag
       |${gcsCopyCommand()} > gsutil_output.txt 2>&1
       |RC_GSUTIL=$$?
       |
       |${recoverWithRequesterPays()}
       |
       |if [ "$$RC_GSUTIL" != "0" ]; then
       |  echo "Failed to download the file. Error: $$(cat gsutil_output.txt)" >&2
       |  exit "$$RC_GSUTIL"
       |else
       |  echo "Download complete!"
       |  exit 0
       |fi
       |""".stripMargin
  }
}

object GcsUriDownloader {
  private final val RequesterPaysErrorMsg = "requester pays bucket but no user project"
}
