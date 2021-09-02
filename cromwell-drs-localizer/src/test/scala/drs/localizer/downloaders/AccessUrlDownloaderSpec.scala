package drs.localizer.downloaders

import cats.effect.ExitCode
import cloud.nio.impl.drs.AccessUrl
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks._
import org.scalatest.prop.TableFor3

class AccessUrlDownloaderSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  it should "return the correct download script for a url-only access URL, no requester pays" in {
    val fakeDownloadLocation = "/root/foo/foo-123.bam"
    val fakeAccessUrl = "http://abc/def/ghi.bam"

    val downloader = AccessUrlDownloader(
      accessUrl = AccessUrl(url = fakeAccessUrl, headers = None),
      downloadLoc = fakeDownloadLocation,
      hashes = None
    )

    val expected = s"""mkdir -p $$(dirname '$fakeDownloadLocation') && rm -f '$fakeDownloadLocation' && getm --checksum-algorithm 'null' --checksum 'null' --filepath '$fakeDownloadLocation' '$fakeAccessUrl'"""

    downloader.generateDownloadScript() shouldBe expected
  }

  val results: TableFor3[Int, String, DownloadResult] = Table(
    ("exitCode", "stderr", "download result"),
    (0, "", DownloadSuccess),
    (0, "what the", RetryableDownloadFailure(ExitCode(0))),
    (0, "oh me oh my: AssertionError: Checksum failed!!!", ChecksumFailure),
    (1, "oh me oh my: AssertionError: Checksum failed!!!", RetryableDownloadFailure(ExitCode(1))),
    (1, " foobar ", RetryableDownloadFailure(ExitCode(1))),
    (0, """ERROR:getm.cli possibly some words "status_code": 503 words""", RetryableDownloadFailure(ExitCode(0))),
    (1, """ERROR:getm.cli possibly some words "status_code": 503 words""", RetryableDownloadFailure(ExitCode(1))),
    (1, """ERROR:getm.cli possibly some words "status_code": 408 more words""", RetryableDownloadFailure(ExitCode(1))),
    (1, """ERROR:getm.cli possibly some words "status_code": 404 even more words""", NonRetryableDownloadFailure(ExitCode(1))),
    (0, """ERROR:getm.cli possibly some words "status_code": 404 even more words""", RetryableDownloadFailure(ExitCode(0))),
  )

  val accessUrlDownloader: AccessUrlDownloader = AccessUrlDownloader(null, null, null)

  forAll(results) { (exitCode, stderr, expected) =>
    it should s"produce $expected for exitCode $exitCode and stderr '$stderr'" in {
      accessUrlDownloader.result(GetmResult(exitCode, null, stderr)) shouldBe expected
    }
  }
}
