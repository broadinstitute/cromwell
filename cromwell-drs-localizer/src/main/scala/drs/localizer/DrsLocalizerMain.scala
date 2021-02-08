package drs.localizer

import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Path}

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO, IOApp}
import cloud.nio.impl.drs.{DrsConfig, MarthaField}
import com.typesafe.scalalogging.StrictLogging

import scala.sys.process._

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

  private final val RequesterPaysErrorMsg = "Bucket is requester pays bucket but no user project provided."
  private final val ExtractGcsUrlErrorMsg = "No resolved url starting with 'gs://' found from Martha response!"

  def getGcsDrsPathResolver: IO[GcsLocalizerDrsPathResolver] = {
    IO {
      val drsConfig = DrsConfig.fromEnv(sys.env)
      new GcsLocalizerDrsPathResolver(drsConfig)
    }
  }


  /*
     Bash to download the GCS file using `gsutil`
   */
  def downloadScript(gcsUrl: String, saJsonPathOption: Option[Path]): String = {

    def gcsCopyCommand(flag: String = ""): String = s"gsutil $flag cp $gcsUrl $downloadLoc"

    def setServiceAccount(): String = {
      saJsonPathOption match {
        case Some(saJsonPath) =>
          s"""# Set gsutil to use the service account returned from Martha
             |gcloud auth activate-service-account --key-file=$saJsonPath > gcloud_output.txt 2>&1
             |RC_GCLOUD=$$?
             |if [ "$$RC_GCLOUD" != "0" ]; then
             |  echo "Failed to activate service account returned from Martha. File won't be downloaded. Error: $$(cat gcloud_output.txt)" >&2
             |  exit "$$RC_GCLOUD"
             |else
             |  echo "Successfully activated service account; Will continue with download. $$(cat gcloud_output.txt)"
             |fi
             |""".stripMargin
        case None => ""
      }
    }

    def recoverWithRequesterPays(): String = {
      requesterPaysProjectIdOption match {
        case Some(userProject) =>
          s"""if [ "$$RC_GSUTIL" != "0" ]; then
             |  # Check if error is requester pays. If yes, retry gsutil copy using project flag
             |  if grep -q '$RequesterPaysErrorMsg' gsutil_output.txt; then
             |    echo "Received 'Bucket is requester pays' error. Attempting again using Requester Pays billing project"
             |    ${gcsCopyCommand(s"-u $userProject")} > gsutil_output.txt 2>&1
             |    RC_GSUTIL=$$?
             |  fi
             |fi""".stripMargin
        case None => ""
      }
    }

    // bash to download the GCS file using gsutil
    s"""set -euo pipefail
       |set +e
       |
       |${setServiceAccount()}
       |
       |# Run gsutil copy without using project flag
       |${gcsCopyCommand()} > gsutil_output.txt 2>&1
       |RC_GSUTIL=$$?
       |
       |${recoverWithRequesterPays()}
       |
       |if [ "$$RC_GSUTIL" != "0" ]; then
       |  echo "Failed to download the file. Error: $$(cat gsutil_output.txt)" >&2
       |  exit "$$RC_GSUTIL"
       |else
       |  echo "Download complete!"
       |  exit 0
       |fi""".stripMargin
  }


  def downloadFileFromGcs(gcsUrl: String, serviceAccountJsonOption: Option[String]) : IO[ExitCode] = {

    logger.info(s"Requester Pays project ID is $requesterPaysProjectIdOption")
    logger.info(s"Attempting to download $gcsUrl to $downloadLoc")

    val copyProcess = serviceAccountJsonOption match {
      case Some(sa) =>
        // if Martha returned a SA, use that SA for gsutil instead of default credentials
        val tempCredentialDir: Path = Files.createTempDirectory("gcloudTemp_").toAbsolutePath
        val saJsonPath: Path = tempCredentialDir.resolve("sa.json")
        Files.write(saJsonPath, sa.getBytes(StandardCharsets.UTF_8))
        val extraEnv = Map("CLOUDSDK_CONFIG" -> tempCredentialDir.toString)
        val copyCommand = Seq("bash", "-c", downloadScript(gcsUrl, Option(saJsonPath)))
        Process(copyCommand, None, extraEnv.toSeq: _*)
      case None =>
        // No SA returned from Martha. gsutil will use the application default credentials.
        val copyCommand = Seq("bash", "-c", downloadScript(gcsUrl, None))
        Process(copyCommand)
    }

    // run the multiple bash script to download file and log stream sent to stdout and stderr using ProcessLogger
    val returnCode = copyProcess ! ProcessLogger(logger.underlying.info, logger.underlying.error)

    IO(ExitCode(returnCode))
  }


  def resolveAndDownload(): IO[ExitCode] = {
    val fields = NonEmptyList.of(MarthaField.GsUri, MarthaField.GoogleServiceAccount)
    for {
      localizerGcsDrsPathResolver <- getGcsDrsPathResolver
      marthaResponse <- localizerGcsDrsPathResolver.resolveDrsThroughMartha(drsUrl, fields)
      // Currently Martha only supports resolving DRS paths to GCS paths
      gcsUrl <- IO.fromEither(marthaResponse.gsUri.toRight(new RuntimeException(ExtractGcsUrlErrorMsg)))
      serviceAccountJsonOption = marthaResponse.googleServiceAccount.map(_.data.spaces2)
      exitState <- downloadFileFromGcs(gcsUrl, serviceAccountJsonOption)
    } yield exitState
  }
}
