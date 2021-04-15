package drs.localizer

import java.nio.file.{Files, Path}
import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.{AccessUrl, DrsConfig, MarthaField, MarthaResponse}
import common.assertion.CromwellTimeoutSpec
import drs.localizer.downloaders.GcsUriDownloader
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


class DrsLocalizerMainSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  val fakeDownloadLocation = "/root/foo/foo-123.bam"
  val fakeRequesterPaysId = "fake-billing-project"

  behavior of "DrsLocalizerMain"

  it should "fail if drs input is not passed" in {
    DrsLocalizerMain.run(List(fakeDownloadLocation)).unsafeRunSync() shouldBe ExitCode.Error
  }

  it should "fail if download location is not passed" in {
    DrsLocalizerMain.run(List(MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly)).unsafeRunSync() shouldBe ExitCode.Error
  }

  it should "accept arguments and run successfully without Requester Pays ID" in {
    val mockDrsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly, fakeDownloadLocation, None)
    val expected = GcsUriDownloader(
      gcsUrl = "gs://abc/foo-123/abc123",
      serviceAccountJson = None,
      downloadLoc = fakeDownloadLocation,
      requesterPaysProjectIdOption = None)
    mockDrsLocalizer.resolve().unsafeRunSync() shouldBe expected
  }

  it should "run successfully with all 3 arguments" in {
    val mockDrsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly, fakeDownloadLocation, Option(fakeRequesterPaysId))
    val expected = GcsUriDownloader(
      gcsUrl = "gs://abc/foo-123/abc123",
      serviceAccountJson = None,
      downloadLoc = fakeDownloadLocation,
      requesterPaysProjectIdOption = Option(fakeRequesterPaysId))
    mockDrsLocalizer.resolve().unsafeRunSync() shouldBe expected
  }

  it should "fail and throw error if Martha response does not have gs:// url" in {
    val mockDrsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithoutAnyResolution, fakeDownloadLocation, None)

    the[RuntimeException] thrownBy {
      mockDrsLocalizer.resolve().unsafeRunSync()
    } should have message "No access URL nor GCS URI starting with 'gs://' found in Martha response!"
  }

  it should "return correct download script for a drs url without Requester Pays ID and Google SA returned from Martha" in {
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
        |fi""".stripMargin

    downloader.gcsDownloadScript(gcsUrl = gcsUrl, saJsonPathOption = None) shouldBe expectedDownloadScript
  }

  it should "inject Requester Pays flag & gcloud auth using SA returned from Martha" in {
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
        |# Set gsutil to use the service account returned from Martha
        |gcloud auth activate-service-account --key-file=${fakeSAJsonPath.toString} > gcloud_output.txt 2>&1
        |RC_GCLOUD=$$?
        |if [ "$$RC_GCLOUD" != "0" ]; then
        |  echo "Failed to activate service account returned from Martha. File won't be downloaded. Error: $$(cat gcloud_output.txt)" >&2
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
        |  if grep -q 'Bucket is requester pays bucket but no user project provided.' gsutil_output.txt; then
        |    echo "Received 'Bucket is requester pays' error. Attempting again using Requester Pays billing project"
        |    gsutil -u fake-billing-project cp $gcsUrl $fakeDownloadLocation > gsutil_output.txt 2>&1
        |    RC_GSUTIL=$$?
        |  fi
        |fi
        |
        |if [ "$$RC_GSUTIL" != "0" ]; then
        |  echo "Failed to download the file. Error: $$(cat gsutil_output.txt)" >&2
        |  exit "$$RC_GSUTIL"
        |else
        |  echo "Download complete!"
        |  exit 0
        |fi""".stripMargin

    downloader.gcsDownloadScript(gcsUrl, Option(fakeSAJsonPath)) shouldBe expectedDownloadScript
  }

  it should "call the correct function for an access url" in {

  }

  it should "return the correct download script for a drs url that resolves to an access URL only" ignore {
    fail("need to write the download script for this")
  }

  it should "return the correct download script for a drs url that resolves to an access URL and a GCS path" ignore {
    fail("need to write the download script for this")
  }

}

object MockDrsPaths {
  val fakeDrsUrlWithGcsResolutionOnly = "drs://abc/foo-123/abc123"
  val fakeDrsUrlWithAccessUrlResolutionOnly = "drs://def/bar-456/def456"
  val fakeDrsUrlWithAccessUrlAndGcsResolution = "drs://ghi/baz-789/ghi789"
  val fakeDrsUrlWithoutAnyResolution = "drs://foo/bar/no-gcs-path"
}


class MockDrsLocalizerMain(drsUrl: String,
                           downloadLoc: String,
                           requesterPaysProjectIdOption: Option[String],
                          )
  extends DrsLocalizerMain(drsUrl, downloadLoc, requesterPaysProjectIdOption) {

  override def getDrsPathResolver: IO[LocalizerDrsPathResolver] = {
    IO {
      new MockLocalizerDrsPathResolver(cloud.nio.impl.drs.MockDrsPaths.mockDrsConfig)
    }
  }
}


class MockLocalizerDrsPathResolver(drsConfig: DrsConfig) extends
  LocalizerDrsPathResolver(drsConfig) {

  override def resolveDrsThroughMartha(drsPath: String, fields: NonEmptyList[MarthaField.Value]): IO[MarthaResponse] = {
    val baseResponse = MarthaResponse(
      size = Option(1234),
      hashes = Option(Map("md5" -> "abc123", "crc32c" -> "34fd67"))
    )

    IO.pure(drsPath) map {
      case MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly =>
        baseResponse.copy(
          gsUri = Option("gs://abc/foo-123/abc123"))
      case MockDrsPaths.fakeDrsUrlWithoutAnyResolution =>
        baseResponse
      case MockDrsPaths.fakeDrsUrlWithAccessUrlResolutionOnly =>
        baseResponse.copy(
          accessUrl = Option(AccessUrl(url = "http://abc/def/ghi.bam", headers = None)))
      case MockDrsPaths.fakeDrsUrlWithAccessUrlAndGcsResolution =>
        baseResponse.copy(
          accessUrl = Option(AccessUrl(url = "http://abc/def/ghi.bam", headers = None)),
          gsUri = Option("gs://some/uri"))
      case _ => throw new RuntimeException("fudge")
    }
  }
}
