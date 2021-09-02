package drs.localizer

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO, IOApp}
import cloud.nio.impl.drs.{AccessUrl, DrsConfig, DrsPathResolver, MarthaField}
import cloud.nio.spi.{CloudNioBackoff, CloudNioSimpleExponentialBackoff}
import com.typesafe.scalalogging.StrictLogging
import drs.localizer.downloaders.AccessUrlDownloader.Hashes
import drs.localizer.downloaders._

import scala.concurrent.duration._
import scala.language.postfixOps

object DrsLocalizerMain extends IOApp with StrictLogging {

  /* This assumes the args are as follows:
      0: DRS input
      1: download location
      2: Optional parameter- Requester Pays Billing project ID
     Martha URL is passed as an environment variable
   */
  override def run(args: List[String]): IO[ExitCode] = {
    val argsLength = args.length

    argsLength match {
      case 2 =>
        new DrsLocalizerMain(args.head, args(1), None).
          resolveAndDownloadWithRetries(downloadRetries = 3, checksumRetries = 1, defaultDownloaderFactory, Option(defaultBackoff)).map(_.exitCode)
      case 3 =>
        new DrsLocalizerMain(args.head, args(1), Option(args(2))).
          resolveAndDownloadWithRetries(downloadRetries = 3, checksumRetries = 1, defaultDownloaderFactory, Option(defaultBackoff)).map(_.exitCode)
      case _ =>
        val argsList = if (args.nonEmpty) args.mkString(",") else "None"
        logger.error(s"Received $argsLength arguments. DRS input and download location path is required. Requester Pays billing project ID is optional. " +
          s"Arguments received: $argsList")
        IO(ExitCode.Error)
    }
  }

  val defaultBackoff: CloudNioBackoff = CloudNioSimpleExponentialBackoff(
    initialInterval = 10 seconds, maxInterval = 60 seconds, multiplier = 2)

  val defaultDownloaderFactory: DownloaderFactory = new DownloaderFactory {
    override def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes): IO[Downloader] =
      IO.pure(AccessUrlDownloader(accessUrl, downloadLoc, hashes))

    override def buildGcsUriDownloader(gcsPath: String, serviceAccountJsonOption: Option[String], downloadLoc: String, requesterPaysProjectOption: Option[String]): IO[Downloader] =
      IO.pure(GcsUriDownloader(gcsPath, serviceAccountJsonOption, downloadLoc, requesterPaysProjectOption))
  }
}

class DrsLocalizerMain(drsUrl: String,
                       downloadLoc: String,
                       requesterPaysProjectIdOption: Option[String]) extends StrictLogging {

  def getDrsPathResolver: IO[DrsLocalizerDrsPathResolver] = {
    IO {
      val drsConfig = DrsConfig.fromEnv(sys.env)
      new DrsLocalizerDrsPathResolver(drsConfig)
    }
  }

  def resolveAndDownloadWithRetries(downloadRetries: Int,
                                    checksumRetries: Int,
                                    downloaderFactory: DownloaderFactory,
                                    backoff: Option[CloudNioBackoff],
                                    downloadAttempt: Int = 0,
                                    checksumAttempt: Int = 0): IO[DownloadResult] = {

    def maybeRetryForChecksumFailure(t: Throwable): IO[DownloadResult] = {
      if (checksumAttempt < checksumRetries) {
        backoff foreach { b => Thread.sleep(b.backoffMillis) }
        logger.warn(s"Attempting retry $checksumAttempt of $checksumRetries checksum retries to download $drsUrl", t)
        resolveAndDownloadWithRetries(downloadRetries, checksumRetries, downloaderFactory, backoff map { _.next }, downloadAttempt, checksumAttempt + 1)
      } else {
        IO.raiseError(new RuntimeException(s"Exhausted $checksumRetries retries to resolve, download and checksum $drsUrl", t))
      }
    }

    def maybeRetryForDownloadFailure(t: Throwable): IO[DownloadResult] = {
      if (downloadAttempt < downloadRetries) {
        backoff foreach { b => Thread.sleep(b.backoffMillis) }
        logger.warn(s"Attempting retry $downloadAttempt of $downloadRetries retries to download $drsUrl", t)
        resolveAndDownloadWithRetries(downloadRetries, checksumRetries, downloaderFactory, backoff map { _.next }, downloadAttempt + 1, checksumAttempt)
      } else {
        IO.raiseError(new RuntimeException(s"Exhausted $downloadRetries retries to resolve, download and checksum $drsUrl", t))
      }
    }

    resolveAndDownload(downloaderFactory).redeemWith({
      maybeRetryForDownloadFailure
    },
    {
      case r: RetryableDownloadFailure =>
        maybeRetryForDownloadFailure(
          new RuntimeException(s"Retryable download error: exit code ${r.exitCode.code} for $drsUrl on retry attempt $downloadAttempt of $downloadRetries"))
      case ChecksumFailure =>
        maybeRetryForChecksumFailure(new RuntimeException(s"Checksum failure for $drsUrl on checksum retry attempt $checksumAttempt of $checksumRetries"))
      case o => IO.pure(o)
    })
  }

  private [localizer] def resolveAndDownload(downloaderFactory: DownloaderFactory): IO[DownloadResult] = {
    resolve(downloaderFactory) flatMap { _.download }
  }

  private [localizer] def resolve(downloaderFactory: DownloaderFactory): IO[Downloader] = {
    val fields = NonEmptyList.of(MarthaField.GsUri, MarthaField.GoogleServiceAccount, MarthaField.AccessUrl, MarthaField.Hashes)
    for {
      resolver <- getDrsPathResolver
      marthaResponse <- resolver.resolveDrsThroughMartha(drsUrl, fields)

      // Currently Martha only supports resolving DRS paths to access URLs or GCS paths.
      downloader <- (marthaResponse.accessUrl, marthaResponse.gsUri) match {
        case (Some(accessUrl), _) =>
          downloaderFactory.buildAccessUrlDownloader(accessUrl, downloadLoc, marthaResponse.hashes)
        case (_, Some(gcsPath)) =>
          val serviceAccountJsonOption = marthaResponse.googleServiceAccount.map(_.data.spaces2)
          downloaderFactory.buildGcsUriDownloader(
            gcsPath = gcsPath,
            serviceAccountJsonOption = serviceAccountJsonOption,
            downloadLoc = downloadLoc,
            requesterPaysProjectOption = requesterPaysProjectIdOption)
        case _ =>
          IO.raiseError(new RuntimeException(DrsPathResolver.ExtractUriErrorMsg))
      }
    } yield downloader
  }
}
