package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.AccessUrl
import com.typesafe.scalalogging.StrictLogging
import common.util.StringUtil._
import drs.localizer.downloaders.AccessUrlDownloader._

import scala.sys.process.{Process, ProcessLogger}
import scala.util.matching.Regex

case class GetmResult(returnCode: Int, stdout: String, stderr: String)

case class AccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes) extends Downloader with StrictLogging {
  def generateDownloadScript(): String = {
    val signedUrl = accessUrl.url
    // TODO headers
    val checksumArgs = GetmChecksum(hashes, accessUrl).args
    s"""mkdir -p $$(dirname '$downloadLoc') && rm -f '$downloadLoc' && getm $checksumArgs --filepath '$downloadLoc' '$signedUrl'"""
  }

  def runGetm: IO[GetmResult] = IO {
    val copyCommand = Seq("bash", "-c", generateDownloadScript())
    val copyProcess = Process(copyCommand)
    val stdout = new StringBuilder()
    val stderr = new StringBuilder()

    val outfn: String => Unit = { s => stdout.append(s); () }
    val errfn: String => Unit = { s => stderr.append(s); () }

    val returnCode = copyProcess ! ProcessLogger(outfn, errfn)

    GetmResult(returnCode, stdout.toString().trim, stderr.toString().trim())
  }

  override def download: IO[DownloadResult] = {
    // We don't want to log the unmasked signed URL here. On a PAPI backend this log will end up under the user's
    // workspace bucket, but that bucket may have visibility different than the data referenced by the signed URL.
    val masked = accessUrl.url.maskSensitiveUri
    logger.info(s"Attempting to download data to '$downloadLoc' from access URL '$masked'.")

    runGetm map toDownloadResult
  }

  def toDownloadResult(getmResult: GetmResult): DownloadResult = {
    getmResult match {
      case GetmResult(0, _, stderr) if stderr.isEmpty =>
        DownloadSuccess
      case GetmResult(0, _, stderr) =>
        stderr match {
          case ChecksumFailureMessage() =>
            ChecksumFailure
          case _ =>
            UnrecognizedRetryableDownloadFailure(ExitCode(0))
        }
      case GetmResult(rc, _, stderr) =>
        stderr match {
          case HttpStatusMessage(status) =>
            val intStatus = Integer.parseInt(status)
            if (intStatus / 100 == 5 || intStatus == 408)
              RecognizedRetryableDownloadFailure(ExitCode(rc))
            else
              NonRetryableDownloadFailure(ExitCode(rc))
          case _ =>
            UnrecognizedRetryableDownloadFailure(ExitCode(rc))
        }
    }
  }
}

object AccessUrlDownloader {
  type Hashes = Option[Map[String, String]]

  val ChecksumFailureMessage: Regex = raw""".*AssertionError: Checksum failed!.*""".r
  val HttpStatusMessage: Regex = raw"""ERROR:getm\.cli.*"status_code":\s*(\d+).*""".r
}
