package drs.localizer

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO}
import cats.syntax.validated._
import cloud.nio.impl.drs.DrsPathResolver.FatalRetryDisposition
import cloud.nio.impl.drs.{AccessUrl, DrsConfig, DrsCredentials, DrsResolverField, DrsResolverResponse}
import common.assertion.CromwellTimeoutSpec
import common.validation.ErrorOr.ErrorOr
import drs.localizer.MockDrsLocalizerDrsPathResolver.{FakeAccessTokenStrategy, FakeHashes}
import drs.localizer.downloaders.AccessUrlDownloader.Hashes
import drs.localizer.downloaders._
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
    mockDrsLocalizer.resolve(DrsLocalizerMain.defaultDownloaderFactory).unsafeRunSync() shouldBe expected
  }

  it should "run successfully with all 3 arguments" in {
    val mockDrsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly, fakeDownloadLocation, Option(fakeRequesterPaysId))
    val expected = GcsUriDownloader(
      gcsUrl = "gs://abc/foo-123/abc123",
      serviceAccountJson = None,
      downloadLoc = fakeDownloadLocation,
      requesterPaysProjectIdOption = Option(fakeRequesterPaysId))
    mockDrsLocalizer.resolve(DrsLocalizerMain.defaultDownloaderFactory).unsafeRunSync() shouldBe expected
  }

  it should "fail and throw error if the DRS Resolver response does not have gs:// url" in {
    val mockDrsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithoutAnyResolution, fakeDownloadLocation, None)

    the[RuntimeException] thrownBy {
      mockDrsLocalizer.resolve(DrsLocalizerMain.defaultDownloaderFactory).unsafeRunSync()
    } should have message "No access URL nor GCS URI starting with 'gs://' found in the DRS Resolver response!"
  }

  it should "resolve to use the correct downloader for an access url" in {
    val mockDrsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithAccessUrlResolutionOnly, fakeDownloadLocation, None)
    val expected = AccessUrlDownloader(
      accessUrl = AccessUrl(url = "http://abc/def/ghi.bam", headers = None),
      downloadLoc = fakeDownloadLocation,
      hashes = FakeHashes
    )
    mockDrsLocalizer.resolve(DrsLocalizerMain.defaultDownloaderFactory).unsafeRunSync() shouldBe expected
  }

  it should "resolve to use the correct downloader for an access url when the DRS Resolver response also contains a gs url" in {
    val mockDrsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithAccessUrlAndGcsResolution, fakeDownloadLocation, None)
    val expected = AccessUrlDownloader(
      accessUrl = AccessUrl(url = "http://abc/def/ghi.bam", headers = None), downloadLoc = fakeDownloadLocation,
      hashes = FakeHashes
    )
    mockDrsLocalizer.resolve(DrsLocalizerMain.defaultDownloaderFactory).unsafeRunSync() shouldBe expected
  }

  it should "not retry on access URL download success" in {
    var actualAttempts = 0

    val drsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithAccessUrlResolutionOnly, fakeDownloadLocation, None) {
      override def resolveAndDownload(downloaderFactory: DownloaderFactory): IO[DownloadResult] = {
        actualAttempts = actualAttempts + 1
        super.resolveAndDownload(downloaderFactory)
      }
    }
    val accessUrlDownloader = IO.pure(new Downloader {
      override def download: IO[DownloadResult] =
        IO.pure(DownloadSuccess)
    })

    val downloaderFactory = new DownloaderFactory {
      override def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes): IO[Downloader] = {
        accessUrlDownloader
      }

      override def buildGcsUriDownloader(gcsPath: String, serviceAccountJsonOption: Option[String], downloadLoc: String, requesterPaysProjectOption: Option[String]): IO[Downloader] = {
        // This test path should never ask for the GCS downloader
        throw new RuntimeException("test failure")
      }
    }

    drsLocalizer.resolveAndDownloadWithRetries(
      downloadRetries = 3,
      checksumRetries = 1,
      downloaderFactory = downloaderFactory,
      backoff = None
    ).unsafeRunSync() shouldBe DownloadSuccess

    actualAttempts shouldBe 1
  }

  it should "retry an appropriate number of times for regular retryable access URL download failures" in {
    var actualAttempts = 0

    val drsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithAccessUrlResolutionOnly, fakeDownloadLocation, None) {
      override def resolveAndDownload(downloaderFactory: DownloaderFactory): IO[DownloadResult] = {
        actualAttempts = actualAttempts + 1
        super.resolveAndDownload(downloaderFactory)
      }
    }
    val accessUrlDownloader = IO.pure(new Downloader {
      override def download: IO[DownloadResult] =
        IO.pure(RecognizedRetryableDownloadFailure(exitCode = ExitCode(0)))
    })

    val downloaderFactory = new DownloaderFactory {
      override def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes): IO[Downloader] = {
        accessUrlDownloader
      }

      override def buildGcsUriDownloader(gcsPath: String, serviceAccountJsonOption: Option[String], downloadLoc: String, requesterPaysProjectOption: Option[String]): IO[Downloader] = {
        // This test path should never ask for the GCS downloader
        throw new RuntimeException("test failure")
      }
    }

    assertThrows[Throwable] {
      drsLocalizer.resolveAndDownloadWithRetries(
        downloadRetries = 3,
        checksumRetries = 1,
        downloaderFactory = downloaderFactory,
        backoff = None
      ).unsafeRunSync()
    }

    actualAttempts shouldBe 4 // 1 initial attempt + 3 retries = 4 total attempts
  }

  it should "retry an appropriate number of times for fatal retryable access URL download failures" in {
    var actualAttempts = 0

    val drsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithAccessUrlResolutionOnly, fakeDownloadLocation, None) {
      override def resolveAndDownload(downloaderFactory: DownloaderFactory): IO[DownloadResult] = {
        actualAttempts = actualAttempts + 1
        IO.raiseError(new RuntimeException("testing: fatal error") with FatalRetryDisposition)
      }
    }

    val accessUrlDownloader = IO.pure(new Downloader {
      override def download: IO[DownloadResult] =
        IO.pure(RecognizedRetryableDownloadFailure(exitCode = ExitCode(0)))
    })

    val downloaderFactory = new DownloaderFactory {
      override def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes): IO[Downloader] = {
        accessUrlDownloader
      }

      override def buildGcsUriDownloader(gcsPath: String, serviceAccountJsonOption: Option[String], downloadLoc: String, requesterPaysProjectOption: Option[String]): IO[Downloader] = {
        // This test path should never ask for the GCS downloader
        throw new RuntimeException("test failure")
      }
    }

    assertThrows[Throwable] {
      drsLocalizer.resolveAndDownloadWithRetries(
        downloadRetries = 3,
        checksumRetries = 1,
        downloaderFactory = downloaderFactory,
        backoff = None
      ).unsafeRunSync()
    }

    actualAttempts shouldBe 1 // 1 and done with a fatal exception
  }

  it should "not retry on GCS URI download success" in {
    var actualAttempts = 0
    val drsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly, fakeDownloadLocation, None) {
      override def resolveAndDownload(downloaderFactory: DownloaderFactory): IO[DownloadResult] = {
        actualAttempts = actualAttempts + 1
        super.resolveAndDownload(downloaderFactory)
      }
    }
    val gcsUriDownloader = IO.pure(new Downloader {
      override def download: IO[DownloadResult] =
        IO.pure(DownloadSuccess)
    })

    val downloaderFactory = new DownloaderFactory {
      override def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes): IO[Downloader] = {
        // This test path should never ask for the access URL downloader
        throw new RuntimeException("test failure")
      }

      override def buildGcsUriDownloader(gcsPath: String, serviceAccountJsonOption: Option[String], downloadLoc: String, requesterPaysProjectOption: Option[String]): IO[Downloader] = {
        gcsUriDownloader
      }
    }

    drsLocalizer.resolveAndDownloadWithRetries(
      downloadRetries = 3,
      checksumRetries = 1,
      downloaderFactory = downloaderFactory,
      backoff = None).unsafeRunSync()

    actualAttempts shouldBe 1
  }

  it should "retry an appropriate number of times for retryable GCS URI download failures" in {
    var actualAttempts = 0
    val drsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly, fakeDownloadLocation, None) {
      override def resolveAndDownload(downloaderFactory: DownloaderFactory): IO[DownloadResult] = {
        actualAttempts = actualAttempts + 1
        super.resolveAndDownload(downloaderFactory)
      }
    }
    val gcsUriDownloader = IO.pure(new Downloader {
      override def download: IO[DownloadResult] =
        IO.pure(RecognizedRetryableDownloadFailure(exitCode = ExitCode(1)))
    })

    val downloaderFactory = new DownloaderFactory {
      override def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes): IO[Downloader] = {
        // This test path should never ask for the access URL downloader
        throw new RuntimeException("test failure")
      }

      override def buildGcsUriDownloader(gcsPath: String, serviceAccountJsonOption: Option[String], downloadLoc: String, requesterPaysProjectOption: Option[String]): IO[Downloader] = {
        gcsUriDownloader
      }
    }

    assertThrows[Throwable] {
      drsLocalizer.resolveAndDownloadWithRetries(
        downloadRetries = 3,
        checksumRetries = 1,
        downloaderFactory = downloaderFactory,
        backoff = None).unsafeRunSync()
    }

    actualAttempts shouldBe 4 // 1 initial attempt + 3 retries = 4 total attempts
  }

  it should "retry an appropriate number of times for checksum failures" in {
    var actualAttempts = 0
    val drsLocalizer = new MockDrsLocalizerMain(MockDrsPaths.fakeDrsUrlWithAccessUrlResolutionOnly, fakeDownloadLocation, None) {
      override def resolveAndDownload(downloaderFactory: DownloaderFactory): IO[DownloadResult] = {
        actualAttempts = actualAttempts + 1
        super.resolveAndDownload(downloaderFactory)
      }
    }
    val accessUrlDownloader = IO.pure(new Downloader {
      override def download: IO[DownloadResult] =
        IO.pure(ChecksumFailure)
    })

    val downloaderFactory = new DownloaderFactory {
      override def buildAccessUrlDownloader(accessUrl: AccessUrl, downloadLoc: String, hashes: Hashes): IO[Downloader] = {
        accessUrlDownloader
      }

      override def buildGcsUriDownloader(gcsPath: String, serviceAccountJsonOption: Option[String], downloadLoc: String, requesterPaysProjectOption: Option[String]): IO[Downloader] = {
        // This test path should never ask for the GCS URI downloader.
        throw new RuntimeException("test failure")
      }
    }

    assertThrows[Throwable] {
      drsLocalizer.resolveAndDownloadWithRetries(
        downloadRetries = 3,
        checksumRetries = 1,
        downloaderFactory = downloaderFactory,
        backoff = None).unsafeRunSync()
    }

    actualAttempts shouldBe 2 // 1 initial attempt + 1 retry = 2 total attempts
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
  extends DrsLocalizerMain(drsUrl, downloadLoc, FakeAccessTokenStrategy, requesterPaysProjectIdOption) {

  override def getDrsPathResolver: IO[DrsLocalizerDrsPathResolver] = {
    IO {
      new MockDrsLocalizerDrsPathResolver(cloud.nio.impl.drs.MockDrsPaths.mockDrsConfig)
    }
  }
}


class MockDrsLocalizerDrsPathResolver(drsConfig: DrsConfig) extends
  DrsLocalizerDrsPathResolver(drsConfig, FakeAccessTokenStrategy) {

  override def resolveDrs(drsPath: String, fields: NonEmptyList[DrsResolverField.Value]): IO[DrsResolverResponse] = {
    val drsResolverResponse = DrsResolverResponse(
      size = Option(1234),
      hashes = FakeHashes
    )

    IO.pure(drsPath) map {
      case MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly =>
        drsResolverResponse.copy(
          gsUri = Option("gs://abc/foo-123/abc123"))
      case MockDrsPaths.fakeDrsUrlWithoutAnyResolution =>
        drsResolverResponse
      case MockDrsPaths.fakeDrsUrlWithAccessUrlResolutionOnly =>
        drsResolverResponse.copy(
          accessUrl = Option(AccessUrl(url = "http://abc/def/ghi.bam", headers = None)))
      case MockDrsPaths.fakeDrsUrlWithAccessUrlAndGcsResolution =>
        drsResolverResponse.copy(
          accessUrl = Option(AccessUrl(url = "http://abc/def/ghi.bam", headers = None)),
          gsUri = Option("gs://some/uri"))
      case e => throw new RuntimeException(s"Unexpected exception in DRS localization test code: $e")
    }
  }
}

object MockDrsLocalizerDrsPathResolver {
  val FakeHashes: Option[Map[String, String]] = Option(Map("md5" -> "abc123", "crc32c" -> "34fd67"))
  val FakeAccessTokenStrategy: DrsCredentials = new DrsCredentials {
    override def getAccessToken: ErrorOr[String] = "testing code: do not call me".invalidNel
  }
}
