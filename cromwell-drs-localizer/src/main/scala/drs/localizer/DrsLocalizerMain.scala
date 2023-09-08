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
          val resolvedUris: IO[List[ResolvedDrsUrl]] = unresolvedUris.flatMap(uris => uris.map(toResolve => resolveDrsUrl(resolver, toResolve)).traverse(identity))

          // Resolve all the URLs, sorting them into buckets by their required downloader type
          val resolvedUrisByType: Map[URIType, List[ResolvedDrsUrl]] = resolvedUris.map { resolvedList =>
            Map(
              URIType.GCS -> resolvedList.filter(resolvedUri => resolvedUri.uriType == URIType.GCS),
              URIType.HTTPS -> resolvedList.filter(resolvedUri => resolvedUri.uriType == URIType.HTTPS),
              URIType.UNKNOWN -> resolvedList.filter(resolvedUri => resolvedUri.uriType == URIType.UNKNOWN)
            )
          }.unsafeRunSync()

          val googleDownload : IO[DownloadResult] = resolvedUrisByType.get(URIType.GCS).map(googleUris => bulkDownloadResolvedGoogleUrls(googleUris, commandLineArguments.googleRequesterPaysProject)).getOrElse(IO.pure(DownloadSuccess))
          val bulkAccessDownload : IO[DownloadResult] = resolvedUrisByType.get(URIType.HTTPS).map(resolvedUris => bulkDownloadResolvedAccessUrls(resolvedUris)).getOrElse(IO.pure(DownloadSuccess))
          resolvedUrisByType.get(URIType.UNKNOWN).map(_ => logger.error("Unrecognized URI format returned from resolver"))

          val googleResult = googleDownload.unsafeRunSync()
          val bulkAccessResult = bulkAccessDownload.unsafeRunSync()
          if(googleResult == DownloadSuccess && bulkAccessResult == DownloadSuccess) ExitCode.Success else ExitCode.Error
        }
        case None => {
          val maybeUrl : Option[UnresolvedDrsUrl] = commandLineArguments.drsObject.flatMap(obj => commandLineArguments.containerPath.map(path => UnresolvedDrsUrl(obj,path)))
          maybeUrl.map(unresolved => {
            val resolved: ResolvedDrsUrl = resolveDrsUrl(resolver, unresolved).unsafeRunSync()
            val downloadResult: DownloadResult = if (resolved.uriType == URIType.HTTPS) bulkDownloadResolvedAccessUrls(List(resolved)).unsafeRunSync() else bulkDownloadResolvedGoogleUrls(List(resolved), commandLineArguments.googleRequesterPaysProject).unsafeRunSync()
            if (downloadResult == DownloadSuccess) ExitCode.Success else ExitCode.Error
          }).getOrElse(ExitCode.Error)
          }
      }
      exitCode
    }

  }

  private def bulkDownloadResolvedAccessUrls(resolvedUrls : List[ResolvedDrsUrl]): IO[DownloadResult] = {
    IO[DownloadResult] {
      if(resolvedUrls.isEmpty) return IO(DownloadSuccess)
      val downloader = defaultDownloaderFactory.buildBulkAccessUrlDownloader(resolvedUrls).unsafeRunSync()
      downloader.download.unsafeRunSync()
    }
  }

  def downloadResolvedGoogleUrl(resolvedUrl : ResolvedDrsUrl, googleRequesterPaysProject : Option[String]) : IO[DownloadResult] = {
    val serviceAccountJsonOption = resolvedUrl.drsResponse.googleServiceAccount.map(_.data.spaces2)
    val downloader = defaultDownloaderFactory.buildGcsUriDownloader(
      gcsPath = resolvedUrl.drsResponse.gsUri.getOrElse("MISSING_GOOGLE_URI"),
      serviceAccountJsonOption = serviceAccountJsonOption,
      downloadLoc = resolvedUrl.downloadDestinationPath,
      requesterPaysProjectOption = googleRequesterPaysProject)
    downloader.flatMap(_.download)
  }
  def bulkDownloadResolvedGoogleUrls(resolvedUrls : List[ResolvedDrsUrl], googleRequesterPaysProject : Option[String]) : IO[DownloadResult] = {
    val downloadIos : List[IO[DownloadResult]] = resolvedUrls.map(url => downloadResolvedGoogleUrl(url, googleRequesterPaysProject))
    val downloadSuccesses : List[DownloadResult] = downloadIos.map(io => io.unsafeRunSync()).filter(res => res==DownloadSuccess)
    val anyFailures = downloadIos.length != downloadSuccesses.length
    if(anyFailures) IO.pure(FatalDownloadFailure(exitCode = ExitCode(1))) else IO.pure(DownloadSuccess)
  }

  //NB: This method calls a synchronous HTTP request
  def resolveDrsUrl(resolverObject: DrsLocalizerDrsPathResolver, drsUrlToResolve: UnresolvedDrsUrl): IO[ResolvedDrsUrl] = {
    IO {
      val fields = NonEmptyList.of(DrsResolverField.GsUri, DrsResolverField.GoogleServiceAccount, DrsResolverField.AccessUrl, DrsResolverField.Hashes)
      //Insert retry logic here.
      val drsResponse = resolverObject.resolveDrs(drsUrlToResolve.drsUrl, fields).unsafeRunSync()
      ResolvedDrsUrl(drsResponse, drsUrlToResolve.downloadDestinationPath, toUriType(drsResponse.accessUrl, drsResponse.gsUri))
    }
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

