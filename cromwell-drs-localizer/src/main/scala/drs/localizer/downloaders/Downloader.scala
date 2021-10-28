package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}

sealed trait DownloadResult {
  def exitCode: ExitCode
}
case object DownloadSuccess extends DownloadResult {
  override def exitCode: ExitCode = ExitCode(0)
}
sealed trait DownloadFailure extends DownloadResult
// For a `FatalDownloadFailure` Cromwell should fail the download immediately with no further download attempts.
case class FatalDownloadFailure(override val exitCode: ExitCode) extends DownloadFailure
sealed trait RetryableDownloadFailure extends DownloadResult
// A regular retryable download failure corresponding to a known set of conditions (http status, stderr content etc.).
case class RecognizedRetryableDownloadFailure(override val exitCode: ExitCode) extends RetryableDownloadFailure
// `UnrecognizedRetryableDownloadFailure` instances are created when the result handling code sees a set of
// circumstances that isn't expected, shrugs and falls back to the default behavior of retrying an attempt.
case class UnrecognizedRetryableDownloadFailure(override val exitCode: ExitCode) extends RetryableDownloadFailure
// Checksum failures represent error-free downloads where a checksum on the downloaded data does not match
// the associated checksum. Checksums failures are currently retried at most once, while `RetryableDownloadFailure`s
// are retried up to three times.
case object ChecksumFailure extends DownloadFailure {
  def exitCode: ExitCode = ExitCode(0)
}

trait Downloader {
  def download: IO[DownloadResult]
}
