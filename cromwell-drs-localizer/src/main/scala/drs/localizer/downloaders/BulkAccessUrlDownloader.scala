package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.AccessUrl
import com.typesafe.scalalogging.StrictLogging
import common.validation.ErrorOr.ErrorOr
import drs.localizer.downloaders.AccessUrlDownloader._

import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Path, Paths}
import scala.sys.process.{Process, ProcessLogger}
import scala.util.Try

case class BulkAccessUrlDownloader(urlToDownloadLocation : Map[AccessUrl, String]) extends Downloader with StrictLogging {

  private def generateBulkDownloadScript: Try[String] = {
    val manifestPath = generateJsonManifest(urlToDownloadLocation).get
    val downloadLoc = "getmDL"
    //val checksumArgs = "pleaseCompile"
    generateJsonManifest(urlToDownloadLocation).map(_ =>  s"""mkdir -p $$(dirname '$downloadLoc') && rm -f '$downloadLoc' && getm --manifest '$manifestPath'""")
  }

  private def toJsonString(accessUrl: String, filepath: String, checksum: ErrorOr[String], checksumAlgorithm: String): String = {
    //NB: trailing comma is being removed in generateJsonManifest
    if (checksum.isValid) {
        s"""  {
           |    "url" : "$accessUrl",
           |    "filepath" : "$filepath",
           |    "checksum" : "$checksum",
           |    "checksum-algorithm" : "$checksumAlgorithm"
           |  },
           |""".stripMargin
    } else {
      s"""  {
         |    "url" : "$accessUrl",
         |    "filepath" : "$filepath"
         |  },
         |""".stripMargin
    }
  }

  def generateJsonManifest(accessUrlToDownloadDest: Map[AccessUrl, String]): Try[Path] = {
    //write a json file that looks like:
    // [
    //  {
    //     "url" : "www.123.com",
    //     "filepath" : "path/to/where/123/should/be/downloaded",
    //     "checksum" : "sdfjndsfjkfsdjsdfkjsdf",
    //     "checksum-algorithm" : "md5"
    //  },
    //  {
    //     "url" : "www.567.com"
    //     "filepath" : "path/to/where/567/should/be/downloaded",
    //     "checksum" : "asdasdasfsdfsdfasdsdfasd",
    //     "checksum-algorithm" : "md5"
    //  }
    // ]
    var jsonString: String = "[\n"
    for ((accessUrl -> downloadDestination) <- accessUrlToDownloadDest) {
      val hashes : Option[Map[String,String]] = Option(Map("" -> ""))
      jsonString += toJsonString(accessUrl.url, downloadDestination, GetmChecksum(hashes, accessUrl).value, GetmChecksum(hashes, accessUrl).getmAlgorithm)
    }
    jsonString = jsonString.substring(0, jsonString.lastIndexOf(",")) //remove trailing comma from array elements
    jsonString += "\n]"
    Try(Files.write(Paths.get("getm-manifest.json"), jsonString.getBytes(StandardCharsets.UTF_8)))
  }

  def runGetm: IO[GetmResult] = {
    val script = generateBulkDownloadScript
    val copyCommand : Seq[String] = Seq("bash", "-c", script.get)
    logger.info(script.get)
    val copyProcess = Process(copyCommand)
    val stderr = new StringBuilder()
    val errorCapture: String => Unit = { s => stderr.append(s); () }
    val returnCode = copyProcess ! ProcessLogger(_ => (), errorCapture)
    logger.info(stderr.toString().trim())
    IO(GetmResult(returnCode, stderr.toString().trim()))
  }

  override def download: IO[DownloadResult] = {
    // We don't want to log the unmasked signed URL here. On a PAPI backend this log will end up under the user's
    // workspace bucket, but that bucket may have visibility different than the data referenced by the signed URL.
    logger.info(s"Attempting to download data")

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
