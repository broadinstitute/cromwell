package drs.localizer

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO, IOApp}
import cats.implicits._
import cloud.nio.impl.drs.DrsPathResolver.{FatalRetryDisposition, RegularRetryDisposition}
import cloud.nio.impl.drs._
import cloud.nio.spi.{CloudNioBackoff, CloudNioSimpleExponentialBackoff}
import com.typesafe.scalalogging.StrictLogging
import drs.localizer.CommandLineParser.AccessTokenStrategy.{Azure, Google}
import drs.localizer.downloaders.AccessUrlDownloader.Hashes
import drs.localizer.downloaders._
import org.apache.commons.csv.{CSVFormat, CSVParser}

import java.io.{File, IOException}
import java.nio.charset.Charset
import scala.concurrent.duration._
import scala.jdk.CollectionConverters._
import scala.language.postfixOps
import spray.json._
import DefaultJsonProtocol._
import drs.localizer.DrsLocalizerMain.{downloadInParallel, generateJsonPathForManifest, loadCSVAsMap, writeJsonManifest}
import drs.localizer.URIType.{GCS, URIType}

import java.nio.file.{Files, Paths}
import java.nio.charset.StandardCharsets
import scala.util.{Failure, Success, Try}

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
    override def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes): IO[Downloader] =
      IO.pure(AccessUrlDownloader(accessUrl, downloadLoc, hashes))

    override def buildGcsUriDownloader(gcsPath: String, serviceAccountJsonOption: Option[String], downloadLoc: String, requesterPaysProjectOption: Option[String]): IO[Downloader] =
      IO.pure(GcsUriDownloader(gcsPath, serviceAccountJsonOption, downloadLoc, requesterPaysProjectOption))
  }

  private def printUsage: IO[ExitCode] = {
    System.err.println(CommandLineParser.Usage)
    IO.pure(ExitCode.Error)
  }

  def runLocalizer(commandLineArguments: CommandLineArguments, drsCredentials: DrsCredentials): IO[ExitCode] = {
    commandLineArguments.manifestPath match {
      case Some(manifestPath) =>
        val manifestFile = new File(manifestPath)
        val csvParser = CSVParser.parse(manifestFile, Charset.defaultCharset(), CSVFormat.DEFAULT)
        val exitCodes: IO[List[ExitCode]] = csvParser.asScala.map(record => {
          val drsObject = record.get(0)
          val containerPath = record.get(1)
          localizeFile(commandLineArguments, drsCredentials, drsObject, containerPath)
        }).toList.sequence
        exitCodes.map(_.find(_ != ExitCode.Success).getOrElse(ExitCode.Success))
      case None =>
        val drsObject = commandLineArguments.drsObject.get
        val containerPath = commandLineArguments.containerPath.get
        localizeFile(commandLineArguments, drsCredentials, drsObject, containerPath)
    }
  }

  private def localizeFile(commandLineArguments: CommandLineArguments, drsCredentials: DrsCredentials, drsObject: String, containerPath: String) = {
    new DrsLocalizerMain(drsObject, containerPath, drsCredentials, commandLineArguments.googleRequesterPaysProject).
      resolveAndDownloadWithRetries(downloadRetries = 3, checksumRetries = 1, defaultDownloaderFactory, Option(defaultBackoff)).map(_.exitCode)
  }

  // Overall plan for bulk DRS downloads:
  // 1) Open the CSV provided to Cromwell Drs Localizer as a map of DRS URI -> Destination Filepath
  // 2) Resolve each DRS URI in the map. This will yield two maps, because a drs URI can resolve to either
  // a gs:// URI or an https:// uri, each of which require different downloaders.
  // 3a) Invoke the GcsUriDownloader (and its underlying gsutil) to download all gs:// URIs in parallel.
  // 3b) Invoke the AccessUrlDownloader (and its underlying getm) to download all https:// URIs in parallel.
  //     - getm requires a --manifest argument, which is a .json file that is a map of https:// URI to destination filepath
  // 4) Report success or failure as granularly as is reasonably possible.
  private def loadCSVAsMap(csvManifestPath : String) : Try[Map[String, String]] = {
    for {
      openFile <- Try(new File(csvManifestPath))
      csvParser <- Try(CSVParser.parse(openFile, Charset.defaultCharset(), CSVFormat.DEFAULT))
      map <- Try(csvParser.getRecords.asScala.map(record => (record.get(0), record.get(1))).toMap)
      if !map.isEmpty //An empty map should be considered an error - we're expecting at least one k/v pair in the manifest.
    } yield map
  }

  //Write string to file. Returns the filepath that the string was written to.
  private def writeJsonManifest(destinationPath : String, jsonData : JsValue): Try[String] = {
    Try(Files.write(Paths.get(destinationPath), jsonData.toString().getBytes(StandardCharsets.UTF_8)).toString)
  }
  // We're going to write a JSON version of the CSV manifest we're provided.
  // Save the new manifest to wherever the CSV is saved, but replace .csv with _manifest.json.
  private def generateJsonPathForManifest(csvPath : String) : String = {
    csvPath.substring(0, csvPath.lastIndexOf('.')) + "_manifest.json"
  }

  // The getm tool can read in a JSON manifest and download all files in parallel.
  // Here, we invoke it with the provided manifest file.
  // The manifest file should have 1 or more key value pairs, where each key is a drs URI
  // and each value is a destination path for where that file should be localized to.
  private def downloadInParallel(jsonManifestPath : String) : IO[ExitCode] = {

  }
}

object URIType extends Enumeration {
  type URIType = Value
  val GCS, HTTPS = Value
}

class DrsLocalizerMain(drsUrl: String,
                       downloadLoc: String,
                       drsCredentials: DrsCredentials,
                       requesterPaysProjectIdOption: Option[String]) extends StrictLogging {

  def getDrsPathResolver: IO[DrsLocalizerDrsPathResolver] = {
    IO {
      val drsConfig = DrsConfig.fromEnv(sys.env)
      logger.info(s"Using ${drsConfig.drsResolverUrl} to resolve DRS Objects")
      new DrsLocalizerDrsPathResolver(drsConfig, drsCredentials)
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
    val fields = NonEmptyList.of(DrsResolverField.GsUri, DrsResolverField.GoogleServiceAccount, DrsResolverField.AccessUrl, DrsResolverField.Hashes)
    for {
      resolver <- getDrsPathResolver
      drsResolverResponse <- resolver.resolveDrs(drsUrl, fields)

      // Currently DRS Resolver only supports resolving DRS paths to access URLs or GCS paths.
      downloader <- (drsResolverResponse.accessUrl, drsResolverResponse.gsUri) match {
        case (Some(accessUrl), _) =>
          downloaderFactory.buildAccessUrlDownloader(accessUrl, downloadLoc, drsResolverResponse.hashes)
        case (_, Some(gcsPath)) =>
          val serviceAccountJsonOption = drsResolverResponse.googleServiceAccount.map(_.data.spaces2)
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

  //TODO: use this fn to resolve drs URIs, split the responses into gs:// and https:// groups, and invoke new downloaders for each of those groups.
  private[localizer] def resolveDrsUrl(resolverObject: IO[DrsLocalizerDrsPathResolver], drsUrlToResolve: String): IO[DrsResolverResponse] = {
    val fields = NonEmptyList.of(DrsResolverField.GsUri, DrsResolverField.GoogleServiceAccount, DrsResolverField.AccessUrl, DrsResolverField.Hashes)
    for {
      resolver <- resolverObject
      drsResolverResponse <- resolver.resolveDrs(drsUrlToResolve, fields)
    } yield drsResolverResponse
  }

  private[localizer] def toUriType(drsResponse: DrsResolverResponse) : Try[URIType] = {
    (drsResponse.accessUrl, drsResponse.gsUri) match {
      case (Some(accessUrl), _) =>
        Try(URIType.HTTPS)
      case (_, Some(gcsPath)) =>
        Try(URIType.GCS)
      case _ =>
        Failure(new RuntimeException(DrsPathResolver.ExtractUriErrorMsg))
    }
  }

  private[localizer] def bulkResolve(): Unit = {
    // Load path resolver once from config
    val resolver = getDrsPathResolver

    // resolve all the URIs once.
    val resolvedUris : Try[Map[DrsResolverResponse, String]] = for {
      drsUriToDownloadLocationMap : Map[String, String] <- loadCSVAsMap("path/to/manifest")
      resolvedUriToDownloadLocationMap : Map[DrsResolverResponse, String] = drsUriToDownloadLocationMap.map(value => (resolveDrsUrl(resolver, value._1), value._2))
    } yield resolvedUriToDownloadLocationMap

    // find the GCS ones and download them
      //find the HTTPS ones and download them
    val exitCode: IO[ExitCode] = for {
      drsResponses <- resolvedUris
      httpsUris = drsResponses.filter(drsResponse => toUriType(drsResponse._1) == URIType.HTTPS)
      uriToDownloadLocation = httpsUris.map(drsResponseToDownloadLoc => (drsResponseToDownloadLoc._1.accessUrl.get.toString, drsResponseToDownloadLoc._2))
      jsonString <- uriToDownloadLocation.toJson
      pathOfJsonFile <- writeJsonManifest(generateJsonPathForManifest("path/to/manifest"), jsonString)
      downloadSuccess <- downloadInParallel(pathOfJsonFile)
    } yield downloadSuccess

  }

}

