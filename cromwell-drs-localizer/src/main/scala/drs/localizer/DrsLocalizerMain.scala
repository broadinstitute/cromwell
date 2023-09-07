package drs.localizer

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO, IOApp}
import cats.implicits.toTraverseOps
import cloud.nio.impl.drs._
import cloud.nio.spi.{CloudNioBackoff, CloudNioSimpleExponentialBackoff}
import com.typesafe.scalalogging.StrictLogging
import drs.localizer.CommandLineParser.AccessTokenStrategy.{Azure, Google}
import drs.localizer.downloaders._
import org.apache.commons.csv.{CSVFormat, CSVParser}

import java.io.File
import java.nio.charset.Charset
import scala.concurrent.duration._
import scala.jdk.CollectionConverters._
import scala.language.postfixOps
import drs.localizer.URIType.URIType

import scala.util.{Failure, Try}
case class UnresolvedDrsUrl(drsUrl: String, downloadDestinationPath: String)
case class ResolvedDrsUrlPendingDownload(drsResponse: DrsResolverResponse, downloadDestinationPath: String, uriType: URIType)
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

    override def buildBulkAccessUrlDownloader(urlsToDownload: List[ResolvedDrsUrlPendingDownload]): IO[Downloader] = {
      IO.pure(BulkAccessUrlDownloader(urlsToDownload))
    }
  }

  private def printUsage: IO[ExitCode] = {
    System.err.println(CommandLineParser.Usage)
    IO.pure(ExitCode.Error)
  }

  private def getDrsPathResolver(drsCredentials: DrsCredentials): IO[DrsLocalizerDrsPathResolver] = {
    IO {
      val drsConfig = DrsConfig.fromEnv(sys.env)
      logger.info(s"Using ${drsConfig.drsResolverUrl} to resolve DRS Objects")
      new DrsLocalizerDrsPathResolver(drsConfig, drsCredentials)
    }
  }

  private[localizer] def toUriType(accessUrl: Option[AccessUrl], gsUri : Option[String]): URIType = {
    (accessUrl, gsUri) match {
      case (_, Some(_)) =>
        URIType.GCS
      case (Some(_), _) =>
        URIType.HTTPS
      case (_,_) =>
        URIType.UNKNOWN
    }
  }
  def runLocalizer(commandLineArguments: CommandLineArguments, drsCredentials: DrsCredentials): IO[ExitCode] = {
    IO {
      // Get a single resolver to resolve URLs with
      val resolver = getDrsPathResolver(drsCredentials).unsafeRunSync()

      val exitCode : ExitCode = commandLineArguments.manifestPath match {
        case Some(manifestPath) => {
          val unresolvedUris: IO[List[UnresolvedDrsUrl]] = loadCSVManifest(manifestPath)
          val resolvedUris: IO[List[ResolvedDrsUrlPendingDownload]] = unresolvedUris.flatMap(uris => uris.map(toResolve => resolveDrsUrl(resolver, toResolve)).traverse(identity))

          // Resolve all the URLs, sorting them into buckets by their required downloader type
          val resolvedUrisByType: Map[URIType, List[ResolvedDrsUrlPendingDownload]] = resolvedUris.map { resolvedList =>
            Map(
              URIType.GCS -> resolvedList.filter(resolvedUri => resolvedUri.uriType == URIType.GCS),
              URIType.HTTPS -> resolvedList.filter(resolvedUri => resolvedUri.uriType == URIType.HTTPS),
              URIType.UNKNOWN -> resolvedList.filter(resolvedUri => resolvedUri.uriType == URIType.UNKNOWN)
            )
          }.unsafeRunSync()

          val googleDownload : IO[List[DownloadResult]] = resolvedUrisByType.get(URIType.GCS).map(googleUris => bulkDownloadResolvedGoogleUrls(googleUris, commandLineArguments.googleRequesterPaysProject)).getOrElse(List())
          val bulkAccessDownload : IO[DownloadResult] = resolvedUrisByType.get(URIType.HTTPS).map(resolvedUris => bulkDownloadResolvedAccessUrls(resolvedUris)).getOrElse(IO(DownloadSuccess))
          resolvedUrisByType.get(URIType.UNKNOWN).map(_ => logger.error("Unrecognized URI format returned from resolver"))

          val bulkAccessResult = bulkAccessDownload.unsafeRunSync()
          val googleResults = googleDownload.unsafeRunSync()
          yield if (googleResults.find(res => res != DownloadSuccess).isEmpty && bulkAccessResult.exitCode == ExitCode.Success) ExitCode.Success else ExitCode.Error
        }
        case None =>
          val drsObject = commandLineArguments.drsObject.get
          val containerPath = commandLineArguments.containerPath.get
          val resolvedURI = resolveDrsUrl(resolver, drsObject)
          val map: Map[DrsResolverResponse, String] = Map((resolvedURI.unsafeRunSync() -> containerPath))
          bulkDownloadResolvedAccessUrls(map, drsCredentials)
        //TODO: Also handle the Google download case here.
      }
    }
  }

  private def bulkDownloadResolvedAccessUrls(resolvedUrls : List[ResolvedDrsUrlPendingDownload]): IO[DownloadResult] = {
    IO[DownloadResult] {
      if(resolvedUrls.isEmpty) return IO(DownloadSuccess)
      val downloader = defaultDownloaderFactory.buildBulkAccessUrlDownloader(resolvedUrls).unsafeRunSync()
      downloader.download.unsafeRunSync()
    }
  }

  def bulkDownloadResolvedGoogleUrls(resolvedUrls : List[ResolvedDrsUrlPendingDownload], googleRequesterPaysProject : Option[String]) : IO[List[DownloadResult]] = {
    resolvedUrls.map(url => {
      val serviceAccountJsonOption = url.drsResponse.googleServiceAccount.map(_.data.spaces2)
      val downloader = defaultDownloaderFactory.buildGcsUriDownloader(
        gcsPath = url.drsResponse.gsUri.getOrElse("MISSING_GOOGLE_URI"),
        serviceAccountJsonOption = serviceAccountJsonOption,
        downloadLoc = url.downloadDestinationPath,
        requesterPaysProjectOption = googleRequesterPaysProject)
      downloader.flatMap(_.download)
    }).traverse(identity)
  }

  //NB: This method calls a synchronous HTTP request
  def resolveDrsUrl(resolverObject: DrsLocalizerDrsPathResolver, drsUrlToResolve: UnresolvedDrsUrl): IO[ResolvedDrsUrlPendingDownload] = {
    IO {
      val fields = NonEmptyList.of(DrsResolverField.GsUri, DrsResolverField.GoogleServiceAccount, DrsResolverField.AccessUrl, DrsResolverField.Hashes)
      //Insert retry logic here.
      val drsResponse = resolverObject.resolveDrs(drsUrlToResolve.drsUrl, fields).unsafeRunSync()
      ResolvedDrsUrlPendingDownload(drsResponse, drsUrlToResolve.downloadDestinationPath, toUriType(drsResponse.accessUrl, drsResponse.gsUri))
    }
  }
  // resolve all the URIs once.
  def resolveCSVManifest(csvManifestPath: String, resolver: DrsLocalizerDrsPathResolver) : Try[Map[DrsResolverResponse, String]] = {
    val drsUriToDownloadLocationMap: Try[Map[String, String]] = loadCSVAsMap(csvManifestPath)
    if (drsUriToDownloadLocationMap.isFailure) return Failure[Map[DrsResolverResponse, String]](new RuntimeException("Could not parse CSV manifest"))
    val resolvedUriToDownloadLocationMap: Map[DrsResolverResponse, String] = drsUriToDownloadLocationMap.get.map(value => (resolveDrsUrl(resolver, value._1).unsafeRunSync(), value._2))
    Try(resolvedUriToDownloadLocationMap)
  }

  /**
    * Helper function to read a CSV file as a map from drs URL to requested download destination
    * @param csvManifestPath Path to a CSV file where each row is something like: drs://asdf.ghj, path/to/my/directory
    */
  private def loadCSVManifest(csvManifestPath : String) : IO[List[UnresolvedDrsUrl]] = {
    IO {
      val openFile = new File(csvManifestPath)
      val csvParser = CSVParser.parse(openFile, Charset.defaultCharset(), CSVFormat.DEFAULT)
      val list = csvParser.getRecords.asScala.map(record => UnresolvedDrsUrl(record.get(0), record.get(1))).toList
      list
    }
  }
}

object URIType extends Enumeration {
  type URIType = Value
  val GCS, HTTPS, UNKNOWN = Value
}

class DrsLocalizerMain(drsUrl: String,
                       downloadLoc: String,
                       drsCredentials: DrsCredentials,
                       requesterPaysProjectIdOption: Option[String]) extends StrictLogging {
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

