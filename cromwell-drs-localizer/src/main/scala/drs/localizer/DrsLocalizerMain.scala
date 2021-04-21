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

  def getDrsPathResolver: IO[DrsLocalizerDrsPathResolver] = {
    IO {
      val drsConfig = DrsConfig.fromEnv(sys.env)
      new DrsLocalizerDrsPathResolver(drsConfig)
    }
  }

  def resolveAndDownload(): IO[ExitCode] = resolve flatMap { _.download }

  def resolve(): IO[Downloader] = {
    /*
    Initially do NOT ask for the accessUrl, but do ask for the accessMethodType to know if we should go back and
    attempt to generate an accessUrl.

    As of 2021-04-20 it is not clear if the underlying accessUrl generators will be able to handle the expected load.
    Separating the request for the accessUrl allows this localizer to fall back to the gsUri if the accessUrl
    generation request fails.

    We do NOT want to just retry generic server errors with an uneducated guess that any and all HTTP 5xx responses
    automatically infer that the accessUrl generator is overloaded. If the localizer were to assume that, then ALL
    overloaded DRS providers would be subjected to extra followup requests, including DRS providers that would _never_
    generate an accessUrl.

    A reminder that ailing DRS providers returning HTTP 5xx errors, or even HTTP 429, are already hammered by
    `martha_v3` rapidly retrying due to WA-90.
     */
    val initialRequestFields =
      NonEmptyList.of(MarthaField.GsUri, MarthaField.GoogleServiceAccount, MarthaField.AccessMethodType)
    val accessUrlFields = NonEmptyList.of(MarthaField.AccessUrl)
    for {
      resolver <- getDrsPathResolver
      initialMarthaResponse <- resolver.resolveDrsThroughMartha(drsUrl, initialRequestFields)
      gsUriOption = initialMarthaResponse.gsUri
      serviceAccountJsonOption = initialMarthaResponse.googleServiceAccount.map(_.data.spaces2)
      accessMethodTypeOption = initialMarthaResponse.accessMethodType

      accessUrlOption <- accessMethodTypeOption match {
        case None =>
          // If there is no accessMethodType reported by martha_v3, do not ask for an accessUrl
          IO.pure(None)
        case Some(_) =>
          // Try to get an accessUrl...
          resolver
            .resolveDrsThroughMartha(drsUrl, accessUrlFields)
            .map(_.accessUrl)
            .handleErrorWith(
              throwable =>
                gsUriOption match {
                  case Some(_) =>
                    // Do not worry about the throwable. Use the gsUri.
                    logger.error(
                      "Will use the gsUri due to an error that occurred while requesting the accessUrl: " +
                        throwable.getMessage,
                      throwable,
                    )
                    IO.pure(None)
                  case None =>
                    // There was no gsUri. Return the original throwable.
                    IO.raiseError(throwable)
                }
            )
      }

      // Currently Martha only supports resolving DRS paths to access URLs or GCS paths.
      downloader <- (accessUrlOption, gsUriOption) match {
        case (Some(accessUrl), _) =>
          IO.pure(AccessUrlDownloader(accessUrl, downloadLoc))
        case (_, Some(gcsPath)) =>
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
