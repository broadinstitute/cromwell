package drs.localizer

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO, IOApp}
import cloud.nio.impl.drs.{DrsConfig, MarthaField}
import com.typesafe.scalalogging.StrictLogging
import drs.localizer.downloaders.{AccessUrlDownloader, Downloader, GcsUriDownloader}

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
      case 2 => new DrsLocalizerMain(args.head, args(1), None).resolveAndDownload()
      case 3 => new DrsLocalizerMain(args.head, args(1), Option(args(2))).resolveAndDownload()
      case _ =>
        val argsList = if (args.nonEmpty) args.mkString(",") else "None"
        logger.error(s"Received $argsLength arguments. DRS input and download location path is required. Requester Pays billing project ID is optional. " +
          s"Arguments received: $argsList")
        IO(ExitCode.Error)
    }
  }
}

class DrsLocalizerMain(drsUrl: String,
                       downloadLoc: String,
                       requesterPaysProjectIdOption: Option[String]) extends StrictLogging {

  private final val ExtractUriErrorMsg = "No access URL nor GCS URI starting with 'gs://' found in Martha response!"

  def getDrsPathResolver: IO[LocalizerDrsPathResolver] = {
    IO {
      val drsConfig = DrsConfig.fromEnv(sys.env)
      new LocalizerDrsPathResolver(drsConfig)
    }
  }

  def resolveAndDownload(): IO[ExitCode] = resolve flatMap { _.download }

  def resolve(): IO[Downloader] = {
    val fields = NonEmptyList.of(MarthaField.GsUri, MarthaField.GoogleServiceAccount, MarthaField.AccessUrl)
    for {
      resolver <- getDrsPathResolver
      marthaResponse <- resolver.resolveDrsThroughMartha(drsUrl, fields)

      // Currently Martha only supports resolving DRS paths to access URLs or GCS paths.
      downloader <- (marthaResponse.accessUrl, marthaResponse.gsUri) match {
        case (Some(accessUrl), _) =>
          IO.pure(AccessUrlDownloader(accessUrl, downloadLoc))
        case (_, Some(gcsPath)) =>
          val serviceAccountJsonOption = marthaResponse.googleServiceAccount.map(_.data.spaces2)
          IO.pure(GcsUriDownloader(
            gcsUrl = gcsPath,
            serviceAccountJson = serviceAccountJsonOption,
            downloadLoc = downloadLoc,
            requesterPaysProjectIdOption = requesterPaysProjectIdOption)
          )
        case _ =>
          IO.raiseError(new RuntimeException(ExtractUriErrorMsg))
      }
    } yield downloader
  }
}
