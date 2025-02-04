package drs.localizer.downloaders

import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.{Files, Path}

class GcsUriDownloaderSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  val fakeDownloadLocation = "/root/foo/foo-123.bam"
  val fakeRequesterPaysId = "fake-billing-project"

  it should "return correct download script for a drs url without Requester Pays ID and Google SA returned from the DRS Resolver" in {
    val gcsUrl = "gs://foo/bar.bam"
    val downloader = new GcsUriDownloader(
      gcsUrl = gcsUrl,
      downloadLoc = fakeDownloadLocation,
      serviceAccountJson = None,
      requesterPaysProjectIdOption = None
    )

    val expectedDownloadScript =
      s"""set -euo pipefail
         |set +e
         |
         |
         |
         |# Run gsutil copy without using project flag
         |gsutil  cp $gcsUrl $fakeDownloadLocation > gsutil_output.txt 2>&1
         |RC_GSUTIL=$$?
         |
         |
         |
         |if [ "$$RC_GSUTIL" != "0" ]; then
         |  echo "Failed to download the file. Error: $$(cat gsutil_output.txt)" >&2
         |  exit "$$RC_GSUTIL"
         |else
         |  echo "Download complete!"
         |  exit 0
         |fi
         |""".stripMargin

    downloader.generateDownloadScript(gcsUrl = gcsUrl, saJsonPathOption = None) shouldBe expectedDownloadScript
  }

  it should "inject Requester Pays flag & gcloud auth using SA returned from the DRS Resolver" in {
    val gcsUrl = "gs://foo/bar.bam"
    val downloader = new GcsUriDownloader(
      gcsUrl = gcsUrl,
      downloadLoc = fakeDownloadLocation,
      requesterPaysProjectIdOption = Option(fakeRequesterPaysId),
      serviceAccountJson = None
    )

    val tempCredentialDir: Path = Files.createTempDirectory("gcloudTemp_").toAbsolutePath
    val fakeSAJsonPath: Path = tempCredentialDir.resolve("sa.json")

    val expectedDownloadScript =
      s"""set -euo pipefail
         |set +e
         |
         |# Set gsutil to use the service account returned from the DRS Resolver
         |gcloud auth activate-service-account --key-file=${fakeSAJsonPath.toString} > gcloud_output.txt 2>&1
         |RC_GCLOUD=$$?
         |if [ "$$RC_GCLOUD" != "0" ]; then
         |  echo "Failed to activate service account returned from the DRS Resolver. File won't be downloaded. Error: $$(cat gcloud_output.txt)" >&2
         |  exit "$$RC_GCLOUD"
         |else
         |  echo "Successfully activated service account; Will continue with download. $$(cat gcloud_output.txt)"
         |fi
         |
         |
         |# Run gsutil copy without using project flag
         |gsutil  cp $gcsUrl $fakeDownloadLocation > gsutil_output.txt 2>&1
         |RC_GSUTIL=$$?
         |
         |if [ "$$RC_GSUTIL" != "0" ]; then
         |  # Check if error is requester pays. If yes, retry gsutil copy using project flag
         |  if grep -q 'requester pays bucket but no user project' gsutil_output.txt; then
         |    echo "Received 'Bucket is requester pays' error. Attempting again using Requester Pays billing project"
         |    gsutil -u fake-billing-project cp $gcsUrl $fakeDownloadLocation > gsutil_output.txt 2>&1
         |    RC_GSUTIL=$$?
         |  fi
         |fi
         |
         |
         |if [ "$$RC_GSUTIL" != "0" ]; then
         |  echo "Failed to download the file. Error: $$(cat gsutil_output.txt)" >&2
         |  exit "$$RC_GSUTIL"
         |else
         |  echo "Download complete!"
         |  exit 0
         |fi
         |""".stripMargin

    downloader.generateDownloadScript(gcsUrl, Option(fakeSAJsonPath)) shouldBe expectedDownloadScript
  }
}
