package drs.localizer.downloaders

import cats.effect.IO
import drs.localizer.ResolvedDrsUrlPendingDownload

trait DownloaderFactory {
  def buildBulkAccessUrlDownloader(urlsToDownload: List[ResolvedDrsUrlPendingDownload]) : IO[Downloader]

  def buildGcsUriDownloader(gcsPath: String,
                            serviceAccountJsonOption: Option[String],
                            downloadLoc: String,
                            requesterPaysProjectOption: Option[String]): IO[Downloader]
}
