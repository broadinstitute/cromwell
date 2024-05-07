package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.AccessUrl
import com.typesafe.scalalogging.StrictLogging

import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Path, Paths}
import scala.sys.process.{Process, ProcessLogger}
import scala.util.matching.Regex
import drs.localizer.ResolvedDrsUrl
import spray.json.DefaultJsonProtocol.{listFormat, mapFormat, StringJsonFormat}
import spray.json._

case class GetmResult(returnCode: Int, stderr: String)

/**
  * Getm is a python tool that is used to download resolved DRS uris quickly and in parallel.
  * This class builds a getm-manifest.json file that it uses for input, and builds/executes a shell command
  * to invoke the Getm tool, which is expected to already be installed in the local environment.
  * @param resolvedUrls
  */
case class BulkAccessUrlDownloader(resolvedUrls: List[ResolvedDrsUrl]) extends Downloader with StrictLogging {

  val getmManifestPath: Path = Paths.get("getm-manifest.json")

  /**
    * Write a json manifest to disk that looks like:
    * // [
    * //  {
    * //     "url" : "www.123.com",
    * //     "filepath" : "path/to/where/123/should/be/downloaded",
    * //     "checksum" : "sdfjndsfjkfsdjsdfkjsdf",
    * //     "checksum-algorithm" : "md5"
    * //  },
    * //  {
    * //     "url" : "www.567.com"
    * //     "filepath" : "path/to/where/567/should/be/downloaded",
    * //     "checksum" : "asdasdasfsdfsdfasdsdfasd",
    * //     "checksum-algorithm" : "md5"
    * //  }
    * // ]
    *
    * @param resolvedUrls
    * @return Filepath of a getm-manifest.json that Getm can use to download multiple files in parallel.
    */
  def generateJsonManifest(resolvedUrls: List[ResolvedDrsUrl]): IO[Path] = {
    def resolvedUrlToJsonMap(resolvedUrl: ResolvedDrsUrl): Map[String, String] = {
      val accessUrl: AccessUrl = resolvedUrl.drsResponse.accessUrl.getOrElse(AccessUrl("missing", None))
      resolvedUrl.drsResponse.hashes
        .map { _ =>
          val checksum =
            GetmChecksum(resolvedUrl.drsResponse.hashes, accessUrl).value.getOrElse("error_calculating_checksum")
          val checksumAlgorithm = GetmChecksum(resolvedUrl.drsResponse.hashes, accessUrl).getmAlgorithm
          Map(
            ("url", accessUrl.url),
            ("filepath", resolvedUrl.downloadDestinationPath),
            ("checksum", checksum),
            ("checksum-algorithm", checksumAlgorithm)
          )
        }
        .getOrElse(
          Map(
            ("url", accessUrl.url),
            ("filepath", resolvedUrl.downloadDestinationPath)
          )
        )
    }

    val jsonArray: String = resolvedUrls.map(resolved => resolvedUrlToJsonMap(resolved)).toJson.prettyPrint
    IO(Files.write(getmManifestPath, jsonArray.getBytes(StandardCharsets.UTF_8)))
  }

  def deleteJsonManifest() =
    Files.deleteIfExists(getmManifestPath)

  def generateGetmCommand(pathToMainfestJson: Path): String =
    s"""timeout 24h getm --manifest ${pathToMainfestJson.toString} -vv"""
  def runGetm: IO[GetmResult] =
    generateJsonManifest(resolvedUrls).flatMap { manifestPath =>
      val script = generateGetmCommand(manifestPath)
      val copyCommand: Seq[String] = Seq("bash", "-c", script)
      logger.info(script)
      val copyProcess = Process(copyCommand)
      val stderr = new StringBuilder()
      val errorCapture: String => Unit = { s => stderr.append(s); () }
      val returnCode = copyProcess ! ProcessLogger(_ => (), errorCapture)
      deleteJsonManifest()
      logger.info(stderr.toString().trim())
      IO(GetmResult(returnCode, stderr.toString().trim()))
    }

  override def download: IO[DownloadResult] = {
    // We don't want to log the unmasked signed URL here. On a PAPI backend this log will end up under the user's
    // workspace bucket, but that bucket may have visibility different than the data referenced by the signed URL.
    logger.info(s"Attempting to download data")

    runGetm map toDownloadResult
  }

  def toDownloadResult(getmResult: GetmResult): DownloadResult =
    getmResult match {
      case GetmResult(0, stderr) if stderr.isEmpty =>
        DownloadSuccess
      case GetmResult(0, stderr) =>
        stderr match {
          case BulkAccessUrlDownloader.ChecksumFailureMessage() =>
            ChecksumFailure
          case _ =>
            UnrecognizedRetryableDownloadFailure(ExitCode(0))
        }
      case GetmResult(rc, stderr) =>
        stderr match {
          case BulkAccessUrlDownloader.HttpStatusMessage(status) =>
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

object BulkAccessUrlDownloader {
  type Hashes = Option[Map[String, String]]

  val ChecksumFailureMessage: Regex = raw""".*AssertionError: Checksum failed!.*""".r
  val HttpStatusMessage: Regex = raw"""ERROR:getm\.cli.*"status_code":\s*(\d+).*""".r
}
