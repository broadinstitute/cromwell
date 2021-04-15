package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.AccessUrl
import com.typesafe.scalalogging.StrictLogging

import scala.sys.process.{Process, ProcessLogger}

case class AccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String) extends Downloader with StrictLogging {

  override def download: IO[ExitCode] = {
    val signedUrl = accessUrl.url
      // TODO probably don't want to log the actual signed URL
      logger.info(s"Attempting to download $signedUrl to $downloadLoc")
      // TODO requester pays, possibly refinements to retry strategy
      val script = s"""mkdir -p $$(dirname '$downloadLoc') && curl --location --retry 3 --retry-connrefused --retry-delay 10 --fail --output '$downloadLoc' '$signedUrl'"""
      val copyCommand = Seq("bash", "-c", script)
      val copyProcess = Process(copyCommand)

      val returnCode = copyProcess ! ProcessLogger(logger.underlying.info, logger.underlying.error)

      IO(ExitCode(returnCode))
  }
}
