package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}

sealed trait DownloadResult {
  def exitCode: ExitCode
}
case object DownloadSuccess extends DownloadResult {
  override def exitCode: ExitCode = ExitCode(0)
}
sealed trait DownloadFailure extends DownloadResult
sealed trait RetryableDownloadFailure extends DownloadResult
case class RecognizedRetryableDownloadFailure(override val exitCode: ExitCode) extends RetryableDownloadFailure
// `UnrecognizedRetryableDownloadFailure` instances are created when the result handling code sees a set of
// circumstances that isn't expected, shrugs and falls back to the default behavior of retrying an attempt.
case class UnrecognizedRetryableDownloadFailure(override val exitCode: ExitCode) extends RetryableDownloadFailure
// `NonRetryableDownloadFailure` is a recognized result handling condition that qualifies for a retry attempt
// (assuming the allowed number of retry attempts has not already been met).
case class NonRetryableDownloadFailure(override val exitCode: ExitCode) extends DownloadFailure
case object ChecksumFailure extends DownloadFailure {
  override def exitCode: ExitCode = ExitCode(0)
}

trait Downloader {
  def download: IO[DownloadResult]
}
