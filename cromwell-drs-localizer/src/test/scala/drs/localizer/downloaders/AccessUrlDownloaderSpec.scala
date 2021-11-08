package drs.localizer.downloaders

import cats.effect.ExitCode
import cats.syntax.validated._
import cloud.nio.impl.drs.AccessUrl
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks._

class AccessUrlDownloaderSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  it should "return the correct download script for a url-only access URL, no requester pays" in {
    val fakeDownloadLocation = "/root/foo/foo-123.bam"
    val fakeAccessUrl = "http://abc/def/ghi.bam"

    val downloader = AccessUrlDownloader(
      accessUrl = AccessUrl(url = fakeAccessUrl, headers = None),
      downloadLoc = fakeDownloadLocation,
      hashes = None
    )

    val expected = s"""mkdir -p $$(dirname '$fakeDownloadLocation') && rm -f '$fakeDownloadLocation' && getm --checksum-algorithm 'null' --checksum null --filepath '$fakeDownloadLocation' '$fakeAccessUrl'""".validNel

    downloader.generateDownloadScript shouldBe expected
  }

  {
    val results = Table(
      ("exitCode", "stderr", "download result"),
      (0, "", DownloadSuccess),
      // In `getm` version 0.0.4 checksum failures currently exit 0.
      (0, "oh me oh my: AssertionError: Checksum failed!!!", ChecksumFailure),
      // Unrecognized because of non-zero exit code without an HTTP status, despite what looks like a checksum failure.
      (1, "oh me oh my: AssertionError: Checksum failed!!!", UnrecognizedRetryableDownloadFailure(ExitCode(1))),
      // Unrecognized because of zero exit status with stderr that does not look like a checksum failure.
      (0, "what the", UnrecognizedRetryableDownloadFailure(ExitCode(0))),
      // Unrecognized because of non-zero exit code without an HTTP status.
      (1, " foobar ", UnrecognizedRetryableDownloadFailure(ExitCode(1))),
      // Unrecognized because of zero exit status with stderr that does not look like a checksum failure.
      (0, """ERROR:getm.cli possibly some words "status_code": 503 words""", UnrecognizedRetryableDownloadFailure(ExitCode(0))),
      // Recognized because of non-zero exit status and an HTTP status.
      (1, """ERROR:getm.cli possibly some words "status_code": 503 words""", RecognizedRetryableDownloadFailure(ExitCode(1))),
      // Recognized because of non-zero exit status and an HTTP status.
      (1, """ERROR:getm.cli possibly some words "status_code": 408 more words""", RecognizedRetryableDownloadFailure(ExitCode(1))),
      // Recognized and non-retryable because of non-zero exit status and 404 HTTP status.
      (1, """ERROR:getm.cli possibly some words "status_code": 404 even more words""", FatalDownloadFailure(ExitCode(1))),
      // Unrecognized because of zero exit status and 404 HTTP status.
      (0, """ERROR:getm.cli possibly some words "status_code": 404 even more words""", UnrecognizedRetryableDownloadFailure(ExitCode(0))),
    )

    val accessUrlDownloader = AccessUrlDownloader(null, null, null)

    forAll(results) { (exitCode, stderr, expected) =>
      it should s"produce $expected for exitCode $exitCode and stderr '$stderr'" in {
        accessUrlDownloader.toDownloadResult(GetmResult(exitCode, stderr)) shouldBe expected
      }
    }
  }
}
