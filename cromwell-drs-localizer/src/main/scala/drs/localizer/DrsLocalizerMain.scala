package drs.localizer

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO, IOApp}
import cats.implicits.toTraverseOps
import cloud.nio.impl.drs._
import cloud.nio.spi.{CloudNioBackoff, CloudNioSimpleExponentialBackoff}
import com.typesafe.scalalogging.StrictLogging
import drs.localizer.CommandLineParser.AccessTokenStrategy.{Azure, Google}
import drs.localizer.DrsLocalizerMain.toUriType
import drs.localizer.downloaders._
import org.apache.commons.csv.{CSVFormat, CSVParser}

import java.io.File
import java.nio.charset.Charset
import scala.concurrent.duration._
import scala.jdk.CollectionConverters._
import scala.language.postfixOps
import drs.localizer.URIType.URIType
case class UnresolvedDrsUrl(drsUrl: String, downloadDestinationPath: String)
case class ResolvedDrsUrl(drsResponse: DrsResolverResponse, downloadDestinationPath: String, uriType: URIType)
object DrsLocalizerMain extends IOApp with StrictLogging {

  override def run(args: List[String]): IO[ExitCode] = {
    val parser = buildParser()

    val parsedArgs = parser.parse(args, CommandLineArguments())

    val localize: Option[IO[ExitCode]] = for {
      pa <- parsedArgs
      run <- pa.accessTokenStrategy.collect {
        case Azure => runLocalizer(pa, AzureDrsCredentials(pa.azureIdentityClientId))
        case Google => runLocalizer(pa, GoogleAppDefaultTokenStrategy)
      }
    } yield run

    localize getOrElse printUsage
  }

  def buildParser(): scopt.OptionParser[CommandLineArguments] = new CommandLineParser()

  val defaultBackoff: CloudNioBackoff = CloudNioSimpleExponentialBackoff(
    initialInterval = 10 seconds, maxInterval = 60 seconds, multiplier = 2)

  val defaultDownloaderFactory: DownloaderFactory = new DownloaderFactory {
    override def buildGcsUriDownloader(gcsPath: String, serviceAccountJsonOption: Option[String], downloadLoc: String, requesterPaysProjectOption: Option[String]): IO[Downloader] =
      IO.pure(GcsUriDownloader(gcsPath, serviceAccountJsonOption, downloadLoc, requesterPaysProjectOption))

    override def buildBulkAccessUrlDownloader(urlsToDownload: List[ResolvedDrsUrl]): IO[Downloader] = {
      IO.pure(BulkAccessUrlDownloader(urlsToDownload))
    }
  }

  private def printUsage: IO[ExitCode] = {
    System.err.println(CommandLineParser.Usage)
    IO.pure(ExitCode.Error)
  }

  /**
    * Helper function to read a CSV file as a map from drs URL to requested download destination
    *
    * @param csvManifestPath Path to a CSV file where each row is something like: drs://asdf.ghj, path/to/my/directory
    */
  def loadCSVManifest(csvManifestPath: String): IO[List[UnresolvedDrsUrl]] = {
    IO {
      val openFile = new File(csvManifestPath)
      val csvParser = CSVParser.parse(openFile, Charset.defaultCharset(), CSVFormat.DEFAULT)
      val list = csvParser.getRecords.asScala.map(record => UnresolvedDrsUrl(record.get(0), record.get(1))).toList
      list
    }
  }


  def runLocalizer(commandLineArguments: CommandLineArguments, drsCredentials: DrsCredentials) : IO[ExitCode] = {
    val urlList : IO[List[UnresolvedDrsUrl]] = commandLineArguments.manifestPath match {
      case Some(manifestPath) =>
        loadCSVManifest(manifestPath)
      case None =>
        IO.pure(List(UnresolvedDrsUrl(commandLineArguments.drsObject.get, commandLineArguments.containerPath.get)))
      }
    IO{
      val main = new DrsLocalizerMain(urlList, defaultDownloaderFactory, drsCredentials, commandLineArguments.googleRequesterPaysProject)
      main.resolveAndDownload().unsafeRunSync().exitCode
    }
    }

  def toUriType(accessUrl: Option[AccessUrl], gsUri: Option[String]): URIType = {
    (accessUrl, gsUri) match {
      case (_, Some(_)) =>
        URIType.GCS
      case (Some(_), _) =>
        URIType.HTTPS
      case (_, _) =>
        URIType.UNKNOWN
    }
  }
  }

object URIType extends Enumeration {
  type URIType = Value
  val GCS, HTTPS, UNKNOWN = Value
}

class DrsLocalizerMain(toResolveAndDownload: IO[List[UnresolvedDrsUrl]],
                       downloaderFactory: DownloaderFactory,
                       drsCredentials: DrsCredentials,
                       requesterPaysProjectIdOption: Option[String]) extends StrictLogging {
  def getDrsPathResolver: IO[DrsLocalizerDrsPathResolver] = {
    IO {
      val drsConfig = DrsConfig.fromEnv(sys.env)
      logger.info(s"Using ${drsConfig.drsResolverUrl} to resolve DRS Objects")
      new DrsLocalizerDrsPathResolver(drsConfig, drsCredentials)
    }
  }

  def resolveUrls(unresolvedUrls: IO[List[UnresolvedDrsUrl]]) : IO[List[ResolvedDrsUrl]] = {
    unresolvedUrls.flatMap{unresolvedList =>
      getDrsPathResolver.flatMap{resolver =>
        unresolvedList.map{unresolvedUrl =>
          resolveSingleUrl(resolver, unresolvedUrl)
        }.traverse(identity)
      }
    }
  }

  def resolveSingleUrl(resolverObject: DrsLocalizerDrsPathResolver, drsUrlToResolve: UnresolvedDrsUrl): IO[ResolvedDrsUrl] = {
    IO {
      val fields = NonEmptyList.of(DrsResolverField.GsUri, DrsResolverField.GoogleServiceAccount, DrsResolverField.AccessUrl, DrsResolverField.Hashes)
      //Insert retry logic here.
      val drsResponse = resolverObject.resolveDrs(drsUrlToResolve.drsUrl, fields).unsafeRunSync()
      ResolvedDrsUrl(drsResponse, drsUrlToResolve.downloadDestinationPath, toUriType(drsResponse.accessUrl, drsResponse.gsUri))
    }
  }



  def resolveAndDownload() : IO[DownloadResult] = {
    IO {
      val downloaders: List[Downloader] = buildDownloaders().unsafeRunSync()
      val results: List[DownloadResult] = downloaders.map(downloader => downloader.download.unsafeRunSync())
      results.find(res => res!= DownloadSuccess).getOrElse(DownloadSuccess)
    }
  }
  def buildDownloaders() : IO[List[Downloader]] = {
    resolveUrls(toResolveAndDownload).flatMap { pendingDownloads =>
      val accessUrls = pendingDownloads.filter(url => url.uriType == URIType.HTTPS)
      val googleUrls = pendingDownloads.filter(url => url.uriType == URIType.GCS)
      pendingDownloads.filter(url => url.uriType == URIType.UNKNOWN).map(
        unknown => logger.error(
          s"Url is not an https:// or  gs:// URL: " +
            s"${unknown.drsResponse.accessUrl.getOrElse(unknown.drsResponse.gsUri.getOrElse("missing_url"))}"))
      val bulkDownloader: Option[List[IO[Downloader]]] = if(accessUrls.isEmpty) None else Option(List(buildBulkAccessUrlDownloader(accessUrls)))
      val googleDownloaders: Option[List[IO[Downloader]]] = if(googleUrls.isEmpty) None else  Option(buildGoogleDownloaders(googleUrls))
      val combined: List[IO[Downloader]] = googleDownloaders.map(list => list).getOrElse(List()) ++ bulkDownloader.map(list => list).getOrElse(List())
      combined.traverse(identity)
    }
  }

  def buildGoogleDownloaders(resolvedGoogleUrls: List[ResolvedDrsUrl]) : List[IO[Downloader]] = {
    resolvedGoogleUrls.map{url=>
      downloaderFactory.buildGcsUriDownloader(
        gcsPath = url.drsResponse.gsUri.get,
        serviceAccountJsonOption = url.drsResponse.googleServiceAccount.map(_.data.spaces2),
        downloadLoc = url.downloadDestinationPath,
        requesterPaysProjectOption = requesterPaysProjectIdOption)
    }
  }
  def buildBulkAccessUrlDownloader(resolvedUrls: List[ResolvedDrsUrl]) : IO[Downloader] = {
    downloaderFactory.buildBulkAccessUrlDownloader(resolvedUrls)
  }


  /*
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
*/
    /*
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
*/
      /*
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
    }*


    IO.raiseError(new RuntimeException(s"Exhausted $downloadRetries download retries to resolve, download and checksum $drsUrl", t))
  }
         */
}

