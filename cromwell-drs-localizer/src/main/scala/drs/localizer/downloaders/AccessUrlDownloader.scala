package drs.localizer.downloaders

import cats.data.Validated.{Invalid, Valid}
import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.AccessUrl

import java.nio.file.{Files, Path, Paths}
import java.nio.charset.StandardCharsets
import com.typesafe.scalalogging.StrictLogging
import common.exception.AggregatedMessageException
import common.util.StringUtil._
import common.validation.ErrorOr.ErrorOr
import drs.localizer.downloaders.AccessUrlDownloader._

import scala.sys.process.{Process, ProcessLogger}
import scala.util.Try
import scala.util.matching.Regex

case class GetmResult(returnCode: Int, stderr: String)

case class AccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes) extends Downloader with StrictLogging {
  def generateDownloadScript: ErrorOr[String] = {
    val signedUrl = accessUrl.url
    GetmChecksum(hashes, accessUrl).args map { checksumArgs =>
      s"""mkdir -p $$(dirname '$downloadLoc') && rm -f '$downloadLoc' && getm $checksumArgs --filepath '$downloadLoc' '$signedUrl'"""
    }
  }

  def generateBulkDownloadScript: ErrorOr[String] = {
    val manifestPath = generateJsonManifest()
  }

  def toJsonString(accessUrl: String, filepath : String, checksum : ErrorOr[String], checksumAlgorithm : String): String = {
    //NB: trailing comma is being removed in generateJsonManifest
    if (checksum.isValid) {
      return s"""{
                | "url" : "$accessUrl"
                | "filepath" : "$filepath"
                | "checksum" : "$checksum"
                | "checksum-algorithm" : "$checksumAlgorithm"
                | },
                |""".stripMargin
    }

    s"""{
         | "url" : "$accessUrl"
         | "filepath" : "$filepath"
         | },
         |""".stripMargin
  }
  def generateJsonManifest(accessUrlToDownloadDest : Map[AccessUrl, String]) : Try[Path] = {
    //write a json file that looks like:
    // [
    //  {
    //     "url" : "www.123.com"
    //     "filepath" : "path/to/where/123/should/be/downloaded"
    //     "checksum" : "sdfjndsfjkfsdjsdfkjsdf"
    //     "checksum-algorithm" : "md5"
    //  },
    //  {
    //     "url" : "www.567.com"
    //     "filepath" : "path/to/where/567/should/be/downloaded"
    //     "checksum" : "asdasdasfsdfsdfasdsdfasd"
    //     "checksum-algorithm" : "md5"
    //  }
    // ]
    var jsonString : String = "[\n"
    for ( (accessUrl -> downloadDestination) <- accessUrlToDownloadDest) {
      jsonString += toJsonString(accessUrl.url, downloadDestination, GetmChecksum(hashes, accessUrl).escapedChecksum, GetmChecksum(hashes, accessUrl).getmAlgorithm)
    }
    jsonString = jsonString.substring(0, jsonString.lastIndexOf(",")) //remove trailing comma from array elements
    jsonString += "]"
    Try(Files.write(Paths.get("getm-manifest.json"), jsonString.getBytes(StandardCharsets.UTF_8)))
  }

  def runGetm: IO[GetmResult] = {
    generateDownloadScript match {
      case Invalid(errors) =>
        IO.raiseError(AggregatedMessageException("Error generating access URL download script", errors.toList))
      case Valid(script) => IO {
        val copyCommand = Seq("bash", "-c", script)
        val copyProcess = Process(copyCommand)

        val stderr = new StringBuilder()
        val errorCapture: String => Unit = { s => stderr.append(s); () }

        // As of `getm` version 0.0.4 the contents of stdout do not appear to be interesting (only a progress bar
        // with no option to suppress it), so ignore stdout for now. If stdout becomes interesting in future versions
        // of `getm` it can be captured just like stderr is being captured here.
        val returnCode = copyProcess ! ProcessLogger(_ => (), errorCapture)

        GetmResult(returnCode, stderr.toString().trim())
      }
    }
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
                RecognizedRetryableDownloadFailure(ExitCode(rc))
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
