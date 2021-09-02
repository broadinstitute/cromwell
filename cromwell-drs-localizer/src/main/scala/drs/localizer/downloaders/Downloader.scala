package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}

sealed trait DownloadResult {
  def exitCode: ExitCode
}
case object DownloadSuccess extends DownloadResult {
  override def exitCode: ExitCode = ExitCode(0)
}
sealed trait DownloadFailure extends DownloadResult
case class RetryableDownloadFailure(override val exitCode: ExitCode) extends DownloadFailure
case class NonRetryableDownloadFailure(override val exitCode: ExitCode) extends DownloadFailure
case object ChecksumFailure extends DownloadFailure {
  override def exitCode: ExitCode = ExitCode(0)
}

trait Downloader {
  def download: IO[DownloadResult]
}
