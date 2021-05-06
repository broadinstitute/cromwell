package drs.localizer.downloaders

import cloud.nio.impl.drs.AccessUrl
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class AccessUrlDownloaderSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  it should "return the correct download script for a url-only access URL, no requester pays" in {
    val fakeDownloadLocation = "/root/foo/foo-123.bam"
    val fakeAccessUrl = "http://abc/def/ghi.bam"

    val downloader = AccessUrlDownloader(
      accessUrl = AccessUrl(url = fakeAccessUrl, headers = None),
      downloadLoc = fakeDownloadLocation
    )

    val expected = s"""mkdir -p $$(dirname '$fakeDownloadLocation') && curl --location --retry 3 --retry-connrefused --retry-delay 10 --fail --output '$fakeDownloadLocation' '$fakeAccessUrl'"""

    downloader.generateDownloadScript() shouldBe expected
  }
}
