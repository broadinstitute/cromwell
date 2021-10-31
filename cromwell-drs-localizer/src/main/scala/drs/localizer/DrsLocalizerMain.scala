package drs.localizer

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO, IOApp}
import cloud.nio.impl.drs.DrsPathResolver.{FatalRetryDisposition, RegularRetryDisposition}
import cloud.nio.impl.drs.{AccessUrl, DrsConfig, DrsPathResolver, MarthaField}
import cloud.nio.spi.{CloudNioBackoff, CloudNioSimpleExponentialBackoff}
import com.typesafe.scalalogging.StrictLogging
import drs.localizer.CommandLineParser.AccessTokenStrategy.{Azure, Google}
import drs.localizer.accesstokens.{AccessTokenStrategy, AzureB2CAccessTokenStrategy, GoogleAccessTokenStrategy}
import drs.localizer.downloaders.AccessUrlDownloader.Hashes
import drs.localizer.downloaders._

import scala.concurrent.duration._
import scala.language.postfixOps

object DrsLocalizerMain extends IOApp with StrictLogging {

  override def run(args: List[String]): IO[ExitCode] = {
    val parser = buildParser()

    val parsedArgs = parser.parse(args, CommandLineArguments())

    val localize: Option[IO[ExitCode]] = for {
      pa <- parsedArgs
      run <- pa.accessTokenStrategy.collect {
        case Azure => runLocalizer(pa, AzureB2CAccessTokenStrategy(pa))
        case Google => runLocalizer(pa, GoogleAccessTokenStrategy)
      }
    } yield run

    localize getOrElse printUsage
  }

  def buildParser(): scopt.OptionParser[CommandLineArguments] = new CommandLineParser()

  val defaultBackoff: CloudNioBackoff = CloudNioSimpleExponentialBackoff(
    initialInterval = 10 seconds, maxInterval = 60 seconds, multiplier = 2)

  val defaultDownloaderFactory: DownloaderFactory = new DownloaderFactory {
    override def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes): IO[Downloader] =
      IO.pure(AccessUrlDownloader(accessUrl, downloadLoc, hashes))

    override def buildGcsUriDownloader(gcsPath: String, serviceAccountJsonOption: Option[String], downloadLoc: String, requesterPaysProjectOption: Option[String]): IO[Downloader] =
      IO.pure(GcsUriDownloader(gcsPath, serviceAccountJsonOption, downloadLoc, requesterPaysProjectOption))
  }

  private def printUsage: IO[ExitCode] = {
    System.err.println(CommandLineParser.Usage)
    IO.pure(ExitCode.Error)
  }

  def runLocalizer(commandLineArguments: CommandLineArguments, accessTokenStrategy: AccessTokenStrategy): IO[ExitCode] = {
    val drsObject = commandLineArguments.drsObject.get
    val containerPath = commandLineArguments.containerPath.get
    new DrsLocalizerMain(drsObject, containerPath, accessTokenStrategy, commandLineArguments.googleRequesterPaysProject).
      resolveAndDownloadWithRetries(downloadRetries = 3, checksumRetries = 1, defaultDownloaderFactory, Option(defaultBackoff)).map(_.exitCode)
  }
}

class DrsLocalizerMain(drsUrl: String,
                       downloadLoc: String,
                       accessTokenStrategy: AccessTokenStrategy,
                       requesterPaysProjectIdOption: Option[String]) extends StrictLogging {

  def getDrsPathResolver: IO[DrsLocalizerDrsPathResolver] = {
    IO {
      val drsConfig = DrsConfig.fromEnv(sys.env)
      new DrsLocalizerDrsPathResolver(drsConfig, accessTokenStrategy)
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
        // In the event of a checksum failure reset the download attempt to zero.
        resolveAndDownloadWithRetries(downloadRetries, checksumRetries, downloaderFactory, backoff map { _.next }, 0, checksumAttempt + 1)
      } else {
        IO.raiseError(new RuntimeException(s"Exhausted $checksumRetries checksum retries to resolve, download and checksum $drsUrl", t))
      }
    }

    def maybeRetryForDownloadFailure(t: Throwable): IO[DownloadResult] = {
      t match {
        case _: FatalRetryDisposition =>
          IO.raiseError(t)
        case _ if downloadAttempt < downloadRetries =>
          backoff foreach { b => Thread.sleep(b.backoffMillis) }
          logger.warn(s"Attempting retry $downloadAttempt of $downloadRetries download retries to download $drsUrl", t)
          resolveAndDownloadWithRetries(downloadRetries, checksumRetries, downloaderFactory, backoff map { _.next }, downloadAttempt + 1, checksumAttempt)
        case _ =>
          IO.raiseError(new RuntimeException(s"Exhausted $downloadRetries download retries to resolve, download and checksum $drsUrl", t))
      }
    }

    resolveAndDownload(downloaderFactory).redeemWith({
      maybeRetryForDownloadFailure
    },
    {
      case f: FatalDownloadFailure =>
        IO.raiseError(new RuntimeException(s"Fatal error downloading DRS object: $f"))
      case r: RetryableDownloadFailure =>
        maybeRetryForDownloadFailure(
          new RuntimeException(s"Retryable download error: $r for $drsUrl on retry attempt $downloadAttempt of $downloadRetries") with RegularRetryDisposition)
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
