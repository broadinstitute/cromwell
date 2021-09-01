package drs.localizer.downloaders

import cats.effect.IO
import cloud.nio.impl.drs.AccessUrl

trait DownloaderFactory {
  type Hashes = Option[Map[String, String]]

  def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes): IO[Downloader]

  def buildGcsUriDownloader(gcsPath: String,
                            serviceAccountJsonOption: Option[String],
                            downloadLoc: String,
                            requesterPaysProjectOption: Option[String]): IO[Downloader]
}
