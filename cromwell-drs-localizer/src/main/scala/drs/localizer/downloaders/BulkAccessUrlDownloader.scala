package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.{AccessUrl, DrsResolverResponse}
import com.typesafe.scalalogging.StrictLogging

import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Path, Paths}
import scala.sys.process.{Process, ProcessLogger}
import scala.util.matching.Regex
import drs.localizer.ResolvedDrsUrl
case class GetmResult(returnCode: Int, stderr: String)

/**
  * Getm is a python tool that is used to download resolved DRS uris quickly and in parallel.
  * This class builds a getm-manifest.json file that it uses for input, and builds/executes a shell command
  * to invoke the Getm tool, which is expected to already be installed in the local environment.
  * @param resolvedUrls
  */
case class BulkAccessUrlDownloader(resolvedUrls : List[ResolvedDrsUrl]) extends Downloader with StrictLogging {
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
  def generateJsonManifest(resolvedUrls : List[ResolvedDrsUrl]): IO[Path] = {
    def toJsonString(drsResponse: DrsResolverResponse, destinationFilepath: String): String = {
      //NB: trailing comma is being removed in generateJsonManifest
      val accessUrl: AccessUrl = drsResponse.accessUrl.getOrElse(AccessUrl("missing", None))
      drsResponse.hashes.map(_ => {
        val checksum = GetmChecksum(drsResponse.hashes, accessUrl).value.getOrElse("error_calculating_checksum")
        val checksumAlgorithm = GetmChecksum(drsResponse.hashes, accessUrl).getmAlgorithm
        s"""  {
           |    "url" : "${accessUrl.url}",
           |    "filepath" : "$destinationFilepath",
           |    "checksum" : "$checksum",
           |    "checksum-algorithm" : "$checksumAlgorithm"
           |  },
           |""".stripMargin
      }).getOrElse(
        s"""  {
           |    "url" : "${accessUrl.url}",
           |    "filepath" : "$destinationFilepath"
           |  },
           |""".stripMargin
      )
    }
    IO {
      var jsonString: String = "[\n"
      for (resolvedUrl <- resolvedUrls) {
        jsonString += toJsonString(resolvedUrl.drsResponse, resolvedUrl.downloadDestinationPath)
      }
      if(jsonString.contains(',')) {
        //remove trailing comma from array elements, but don't crash on empty list.
        jsonString = jsonString.substring(0, jsonString.lastIndexOf(","))
      }
      jsonString += "\n]"
      logger.info(jsonString)
      Files.write(Paths.get("getm-manifest.json"), jsonString.getBytes(StandardCharsets.UTF_8))
    }
  }

  def generateGetmCommand(pathToMainfestJson : Path) : String = {
    s"""getm --manifest ${pathToMainfestJson.toString}"""
  }
  def runGetm: IO[GetmResult] = {
    generateJsonManifest(resolvedUrls).flatMap{ manifestPath =>
      //val script = s"""mkdir -p $$(dirname '$downloadLoc') && rm -f '$downloadLoc' && getm --manifest '$manifestPath'""" //TODO: Check if getm will automatically create directories, or if we need to do it for each file.
      // also consider deleting files already there to make retires a little simpler?

      val script = generateGetmCommand(manifestPath)
      val copyCommand : Seq[String] = Seq("bash", "-c", script)
      logger.info(script)
      val copyProcess = Process(copyCommand)
      val stderr = new StringBuilder()
      val errorCapture: String => Unit = { s => stderr.append(s); () }
      val returnCode = copyProcess ! ProcessLogger(_ => (), errorCapture)
      logger.info(stderr.toString().trim())
      IO(GetmResult(returnCode, stderr.toString().trim()))
    }
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
}

object BulkAccessUrlDownloader{
  type Hashes = Option[Map[String, String]]

  val ChecksumFailureMessage: Regex = raw""".*AssertionError: Checksum failed!.*""".r
  val HttpStatusMessage: Regex = raw"""ERROR:getm\.cli.*"status_code":\s*(\d+).*""".r
}
