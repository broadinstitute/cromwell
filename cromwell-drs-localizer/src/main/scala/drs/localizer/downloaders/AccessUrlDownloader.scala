package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.AccessUrl
import com.typesafe.scalalogging.StrictLogging
import common.util.StringUtil._
import drs.localizer.downloaders.AccessUrlDownloader._

import scala.sys.process.{Process, ProcessLogger}
import scala.util.matching.Regex

case class GetmResult(returnCode: Int, stderr: String)

case class AccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes) extends Downloader with StrictLogging {
  def generateDownloadScript(): String = {
    val signedUrl = accessUrl.url
    val checksumArgs = GetmChecksum(hashes, accessUrl).args
    s"""mkdir -p $$(dirname '$downloadLoc') && rm -f '$downloadLoc' && getm $checksumArgs --filepath '$downloadLoc' '$signedUrl'"""
  }

  def runGetm: IO[GetmResult] = IO {
    val copyCommand = Seq("bash", "-c", generateDownloadScript())
    val copyProcess = Process(copyCommand)

    val stderr = new StringBuilder()
    val errorCapture: String => Unit = { s => stderr.append(s); () }

    // As of `getm` version 0.0.4 the contents of stdout do not appear to be interesting (only a progress bar
    // with no option to suppress it), so ignore stdout for now. If stdout becomes interesting in future versions
    // of `getm` it can be captured just like stderr is being captured here.
    val returnCode = copyProcess ! ProcessLogger(_ => (), errorCapture)

    GetmResult(returnCode, stderr.toString().trim())
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
      case GetmResult(0, stderr) if stderr.isEmpty =>
        DownloadSuccess
      case GetmResult(0, stderr) =>
        stderr match {
          case ChecksumFailureMessage() =>
            ChecksumFailure
          case _ =>
            UnrecognizedRetryableDownloadFailure(ExitCode(0))
        }
      case GetmResult(rc, stderr) =>
        stderr match {
          case HttpStatusMessage(status) =>
            Integer.parseInt(status) match {
              case 408 | 429 =>
                TransientRetryableDownloadFailure(ExitCode(rc))
              case s if s / 100 == 4 =>
                FatalDownloadFailure(ExitCode(rc))
              case s if s / 100 == 5 =>
                RecognizedRetryableDownloadFailure(ExitCode(rc))
              case _ =>
                UnrecognizedRetryableDownloadFailure(ExitCode(rc))
            }
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
