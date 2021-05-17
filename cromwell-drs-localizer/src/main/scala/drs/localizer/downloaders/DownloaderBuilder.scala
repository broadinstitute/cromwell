package drs.localizer.downloaders

import cats.effect.IO
import cloud.nio.impl.drs.AccessUrl

trait DownloaderBuilder {
  def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String): IO[Downloader]

  def buildGcsUriDownloader(gcsPath: String,
                            serviceAccountJsonOption: Option[String],
                            downloadLoc: String,
                            requesterPaysProjectOption: Option[String]): IO[Downloader]
}
