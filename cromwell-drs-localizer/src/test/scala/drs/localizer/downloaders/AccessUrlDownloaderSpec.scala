package drs.localizer.downloaders

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

    val expected = s"""mkdir -p $$(dirname '$fakeDownloadLocation') && rm -f '$fakeDownloadLocation' && curl --silent --write-out '%{http_code}' --location --fail --output '$fakeDownloadLocation' '$fakeAccessUrl'"""

    downloader.generateDownloadScript() shouldBe expected
  }

  val results: TableFor3[Int, String, DownloadResult] = Table(
    ("exitCode", "stderr", "expected"),
    (0, "FIXME", null),
    (0, "FIXME", null),
    (1, "FIXME", null),
    (1, "FIXME", null),
    (1, "FIXME", null),
    (1, "FIXME", null),
    (0, "FIXME", null),
  )

  val accessUrlDownloader: AccessUrlDownloader = AccessUrlDownloader(null, null, null)

  forAll(results) { (exitCode, stderr, expected) =>
    it should s"produce $expected for exitCode $exitCode and stderr '$stderr'" in {
      accessUrlDownloader.result(GetmResult(exitCode, "", "")) shouldBe expected
    }
  }
}
