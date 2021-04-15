package drs.localizer

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.{AccessUrl, DrsConfig, MarthaField, MarthaResponse}
import common.assertion.CromwellTimeoutSpec
import drs.localizer.downloaders.{AccessUrlDownloader, GcsUriDownloader}
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

  it should "resolve to use the correct downloader for a gs uri" in {
    val mockDrsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly, fakeDownloadLocation, None)
    val expected = GcsUriDownloader(
      gcsUrl = "gs://abc/foo-123/abc123",
      downloadLoc = fakeDownloadLocation,
      serviceAccountJson = None,
      requesterPaysProjectIdOption = None
    )
    mockDrsLocalizer.resolve().unsafeRunSync() shouldBe expected
  }

  it should "resolve to use the correct downloader for an access url" in {
    val mockDrsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithAccessUrlResolutionOnly, fakeDownloadLocation, None)
    val expected = AccessUrlDownloader(
      accessUrl = AccessUrl(url = "http://abc/def/ghi.bam", headers = None), downloadLoc = fakeDownloadLocation
    )
    mockDrsLocalizer.resolve().unsafeRunSync() shouldBe expected
  }

  it should "resolve to use the correct downloader for an access url when the Martha response also contains a gs url" in {
    val mockDrsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithAccessUrlAndGcsResolution, fakeDownloadLocation, None)
    val expected = AccessUrlDownloader(
      accessUrl = AccessUrl(url = "http://abc/def/ghi.bam", headers = None), downloadLoc = fakeDownloadLocation
    )
    mockDrsLocalizer.resolve().unsafeRunSync() shouldBe expected
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
    val marthaResponse = MarthaResponse(
      size = Option(1234),
      hashes = Option(Map("md5" -> "abc123", "crc32c" -> "34fd67"))
    )

    IO.pure(drsPath) map {
      case MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly =>
        marthaResponse.copy(
          gsUri = Option("gs://abc/foo-123/abc123"))
      case MockDrsPaths.fakeDrsUrlWithoutAnyResolution =>
        marthaResponse
      case MockDrsPaths.fakeDrsUrlWithAccessUrlResolutionOnly =>
        marthaResponse.copy(
          accessUrl = Option(AccessUrl(url = "http://abc/def/ghi.bam", headers = None)))
      case MockDrsPaths.fakeDrsUrlWithAccessUrlAndGcsResolution =>
        marthaResponse.copy(
          accessUrl = Option(AccessUrl(url = "http://abc/def/ghi.bam", headers = None)),
          gsUri = Option("gs://some/uri"))
      case _ => throw new RuntimeException("fudge")
    }
  }
}
