package drs.localizer.downloaders

import drs.localizer.ResolvedDrsUrl

trait DownloaderFactory {
  def buildBulkAccessUrlDownloader(urlsToDownload: List[ResolvedDrsUrl]): Downloader

  def buildGcsUriDownloader(gcsPath: String,
                            serviceAccountJsonOption: Option[String],
                            downloadLoc: String,
                            requesterPaysProjectOption: Option[String]
  ): Downloader
}
