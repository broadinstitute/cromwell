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
      downloadLoc = fakeDownloadLocation
    )

    val expected = s"""mkdir -p $$(dirname '$fakeDownloadLocation') && rm -f '$fakeDownloadLocation' && curl --silent --write-out '%{http_code}' --location --fail --output '$fakeDownloadLocation' '$fakeAccessUrl'"""

    downloader.generateDownloadScript() shouldBe expected
  }

  val results: TableFor3[Int, String, DownloadResult] = Table(
    ("exitCode", "httpStatus", "download result"),
    (0, "200", DownloadSuccess),
    (0, "408", RetryableDownloadFailure(ExitCode(0))),
    (1, "408", RetryableDownloadFailure(ExitCode(1))),
    (1, " foobar ", RetryableDownloadFailure(ExitCode(1))),
    (1, "503", RetryableDownloadFailure(ExitCode(1))),
    (1, "429", NonRetryableDownloadFailure(ExitCode(1))),
    (7, "429", RetryableDownloadFailure(ExitCode(7))),
    (0, "  ", RetryableDownloadFailure(ExitCode(0))),
  )

  val accessUrlDownloader: AccessUrlDownloader = AccessUrlDownloader(null, null)

  forAll(results) { (exitCode, httpStatus, expected) =>
    it should s"produce $expected for exitCode $exitCode and http status '$httpStatus'" in {
      accessUrlDownloader.result(exitCode, httpStatus, "drs://foo/bar") shouldBe expected
    }
  }
}
