package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.AccessUrl
import com.typesafe.scalalogging.StrictLogging
import common.util.StringUtil._
import drs.localizer.DrsLocalizerMain.defaultDownloaderFactory.Hashes

import scala.sys.process.{Process, ProcessLogger}

case class GetmChecksum(algorithm: String, value: String) {
  def args: String = s"--checksum-algorithm '$algorithm' --checksum '$value''"
}

case class AccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes) extends Downloader with StrictLogging {

  // Return a `getm` checksum algorithm and value or None if no checksum should be checked.
  private [downloaders] def checksumAlgorithmAndValue: Option[GetmChecksum] = {
    hashes match {
      case Some(hashes) if hashes.nonEmpty =>
        // `hashes` uses the Martha keys for these hash algorithms, which in turn are forwarded DRS providers' keys for
        // the algorithms. `getm` has its own notions of what these algorithms are called.
        if (hashes.contains("md5")) {
          hashes.get("md5").map(GetmChecksum("md5", _))
        }
        else if (hashes.contains("crc32c")) {
          hashes.get("crc32c").map(GetmChecksum("gs_crc32c", _))
        }
        // etags could be anything; only ask `getm` to check s3 etags if this actually looks like an s3 signed url.
        else if (hashes.contains("etag") && accessUrl.url.matches("^https://[^/]+\\.s3\\.amazonaws.com/.*")) {
          hashes.get("etag") map (GetmChecksum("s3_etag", _))
        }
        // not pictured: sha256, which does appear in Martha test data but is not currently supported by `getm`.
        else {
          // If this code were running in Cromwell this condition would probably merit a warning but the localizer
          // runs on the VM and at best can only complain to stderr. The `getm` algorithm of `null` is specified which
          // means "do not validate checksums" with the stringified contents of the hashes map as a value.
          Option(GetmChecksum("null", hashes.toString()))
        }
      case _ => None
    }
  }

  def generateDownloadScript(): String = {
    val signedUrl = accessUrl.url
    // TODO headers
    // s"""mkdir -p $$(dirname '$downloadLoc') && rm -f '$downloadLoc' && curl --silent --write-out '%{http_code}' --location --fail --output '$downloadLoc' '$signedUrl'"""
    val checksumArgs = checksumAlgorithmAndValue map { _.args } getOrElse ""
    s"""mkdir -p $$(dirname '$downloadLoc') && rm -f '$downloadLoc' && getm --filepath '$downloadLoc' $checksumArgs '$signedUrl'"""
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
