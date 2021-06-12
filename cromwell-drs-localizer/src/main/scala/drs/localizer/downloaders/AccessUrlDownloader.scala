package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.AccessUrl
import com.typesafe.scalalogging.StrictLogging
import common.util.StringUtil._

import scala.sys.process.{Process, ProcessLogger}

case class AccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String) extends Downloader with StrictLogging {

  def generateDownloadScript(): String = {
    val signedUrl = accessUrl.url
    // TODO headers
    s"""mkdir -p $$(dirname '$downloadLoc') && rm -f '$downloadLoc' && curl --silent --write-out '%{http_code}' --location --fail --output '$downloadLoc' '$signedUrl'"""
  }

  def returnCodeAndHttpStatus: IO[(Int, String)] = IO {
    val copyCommand = Seq("bash", "-c", generateDownloadScript())
    val copyProcess = Process(copyCommand)
    val statusString = new StringBuilder()

    val returnCode = copyProcess ! ProcessLogger({ s => statusString.append(s); () }, logger.underlying.error)
    val httpStatus = statusString.toString().trim

    (returnCode, httpStatus)
  }

  override def download: IO[DownloadResult] = {
    // We don't want to log the unmasked signed URL here. On a PAPI backend this log will end up under the user's
    // workspace bucket, but that bucket may have visibility different than the data referenced by the signed URL.
    val masked = accessUrl.url.maskSensitiveUri
    logger.info(s"Attempting to download data to '$downloadLoc' from access URL '$masked'.")

    for {
      rh <- returnCodeAndHttpStatus
      (returnCode, httpStatus) = rh
      ret = result(returnCode, httpStatus, masked)
    } yield ret
  }

  def result(returnCode: Int, httpStatus: String, maskedUrl: String): DownloadResult = {
    (returnCode, httpStatus) match {
      case (i, _) if i == 7 => // ECONNREFUSED
        RetryableDownloadFailure(exitCode = ExitCode(i))
      case (i, s) =>
        try {
          val intStatus = Integer.parseInt(s)
          if (intStatus / 100 == 2) {
            // Not expecting 3xx here since curl is invoked with --location.
            DownloadSuccess
          } else if (intStatus / 100 == 5 || intStatus == 408)
            // The curl --retry behavior is to retry only 5xx or 408.
            RetryableDownloadFailure(ExitCode(i))
          else
            NonRetryableDownloadFailure(ExitCode(i))
        } catch {
          case e: NumberFormatException =>
            logger.warn(s"Unexpected non-numeric http status code for '$maskedUrl': $httpStatus", e)
            RetryableDownloadFailure(ExitCode(i))
        }
    }
  }
}
