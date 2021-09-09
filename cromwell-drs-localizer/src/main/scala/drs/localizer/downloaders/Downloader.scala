package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}

sealed trait DownloadResult {
  def exitCode: ExitCode
}
case object DownloadSuccess extends DownloadResult {
  override def exitCode: ExitCode = ExitCode(0)
}
sealed trait DownloadFailure extends DownloadResult
// For a `FatalDownloadFailure` Cromwell should not make any further download attempts.
case class FatalDownloadFailure(override val exitCode: ExitCode) extends DownloadResult
sealed trait RetryableDownloadFailure extends DownloadResult
// For a `TransientRetryableDownloadFailure` Cromwell should not increment the retry counter.
case class TransientRetryableDownloadFailure(override val exitCode: ExitCode) extends RetryableDownloadFailure
// A regular retryable download failure
case class RecognizedRetryableDownloadFailure(override val exitCode: ExitCode) extends RetryableDownloadFailure
// `UnrecognizedRetryableDownloadFailure` instances are created when the result handling code sees a set of
// circumstances that isn't expected, shrugs and falls back to the default behavior of retrying an attempt.
// In terms of retries this should be handled identically to a `RetryableDownloadFailure`.
case class UnrecognizedRetryableDownloadFailure(override val exitCode: ExitCode) extends RetryableDownloadFailure
// Checksum failures represent error-free downloads where a checksum on the downloaded data does not match
// the associated checksum.
case object ChecksumFailure extends DownloadFailure {
  def exitCode: ExitCode = ExitCode(0)
}

trait Downloader {
  def download: IO[DownloadResult]
}
