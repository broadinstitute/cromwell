package drs.localizer.downloaders

import cats.effect.IO
import drs.localizer.ResolvedDrsUrl

trait DownloaderFactory {
  def buildBulkAccessUrlDownloader(urlsToDownload: List[ResolvedDrsUrl]) : IO[Downloader]

  def buildGcsUriDownloader(gcsPath: String,
                            serviceAccountJsonOption: Option[String],
                            downloadLoc: String,
                            requesterPaysProjectOption: Option[String]): IO[Downloader]
}
