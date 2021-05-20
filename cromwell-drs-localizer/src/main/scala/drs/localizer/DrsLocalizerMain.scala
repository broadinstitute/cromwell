package drs.localizer

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO, IOApp}
import cloud.nio.impl.drs.{AccessUrl, DrsConfig, DrsPathResolver, MarthaField}
import cloud.nio.spi.{CloudNioBackoff, CloudNioSimpleExponentialBackoff}
import com.typesafe.scalalogging.StrictLogging
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
          resolveAndDownloadWithRetries(retries = 3, defaultDownloadBuilder, Option(defaultBackoff)).map(_.exitCode)
      case 3 =>
        new DrsLocalizerMain(args.head, args(1), Option(args(2))).
          resolveAndDownloadWithRetries(retries = 3, defaultDownloadBuilder, Option(defaultBackoff)).map(_.exitCode)
      case _ =>
        val argsList = if (args.nonEmpty) args.mkString(",") else "None"
        logger.error(s"Received $argsLength arguments. DRS input and download location path is required. Requester Pays billing project ID is optional. " +
          s"Arguments received: $argsList")
        IO(ExitCode.Error)
    }
  }

  val defaultBackoff: CloudNioBackoff = CloudNioSimpleExponentialBackoff(
    initialInterval = 10 seconds, maxInterval = 60 seconds, multiplier = 2)

  val defaultDownloadBuilder: DownloaderBuilder = new DownloaderBuilder {
    override def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String): IO[Downloader] =
      IO.pure(AccessUrlDownloader(accessUrl, downloadLoc))

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

  def resolveAndDownloadWithRetries(retries: Int, downloaderBuilder: DownloaderBuilder,
                                    backoff: Option[CloudNioBackoff],
                                    attempt: Int = 0): IO[DownloadResult] = {

    def maybeResolveAndDownloadWithRetry(t: Throwable): IO[DownloadResult] = {
      if (attempt < retries) {
        backoff foreach { b => Thread.sleep(b.backoffMillis) }
        logger.warn(s"Attempting retry $attempt of $retries retries to download $drsUrl", t)
        resolveAndDownloadWithRetries(retries, downloaderBuilder, backoff map { _.next }, attempt + 1)
      } else {
        IO.raiseError(new RuntimeException(s"Exhausted $retries retries to resolve and download $drsUrl", t))
      }
    }

    resolveAndDownload(downloaderBuilder).redeemWith({
      maybeResolveAndDownloadWithRetry
    },
    {
      case r: RetryableDownloadFailure =>
        maybeResolveAndDownloadWithRetry(
          new RuntimeException(s"Retryable download error: exit code ${r.exitCode.code} for $drsUrl on retry attempt $attempt of $retries"))
      case o => IO.pure(o)
    })
  }

  private [localizer] def resolveAndDownload(downloaderBuilder: DownloaderBuilder): IO[DownloadResult] = {
    resolve(downloaderBuilder) flatMap { _.download }
  }

  private [localizer] def resolve(downloaderBuilder: DownloaderBuilder): IO[Downloader] = {
    val fields = NonEmptyList.of(MarthaField.GsUri, MarthaField.GoogleServiceAccount, MarthaField.AccessUrl)
    for {
      resolver <- getDrsPathResolver
      marthaResponse <- resolver.resolveDrsThroughMartha(drsUrl, fields)

      // Currently Martha only supports resolving DRS paths to access URLs or GCS paths.
      downloader <- (marthaResponse.accessUrl, marthaResponse.gsUri) match {
        case (Some(accessUrl), _) =>
          downloaderBuilder.buildAccessUrlDownloader(accessUrl, downloadLoc)
        case (_, Some(gcsPath)) =>
          val serviceAccountJsonOption = marthaResponse.googleServiceAccount.map(_.data.spaces2)
          downloaderBuilder.buildGcsUriDownloader(
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
