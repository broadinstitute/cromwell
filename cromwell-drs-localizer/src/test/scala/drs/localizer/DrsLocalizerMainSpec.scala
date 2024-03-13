package drs.localizer

import cats.data.NonEmptyList
import cats.effect.{ExitCode, IO}
import cats.syntax.validated._
import drs.localizer.MockDrsPaths.{fakeAccessUrls, fakeDrsUrlWithGcsResolutionOnly, fakeGoogleUrls}
import cloud.nio.impl.drs.{AccessUrl, DrsConfig, DrsCredentials, DrsPathResolver, DrsResolverField, DrsResolverResponse}
import common.assertion.CromwellTimeoutSpec
import common.validation.ErrorOr.ErrorOr
import drs.localizer.MockDrsLocalizerDrsPathResolver.{FakeAccessTokenStrategy, FakeHashes}
import drs.localizer.downloaders._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

class DrsLocalizerMainSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  val fakeDownloadLocation = "/root/foo/foo-123.bam"
  val fakeRequesterPaysId = "fake-billing-project"

  val fakeGoogleInput: IO[List[UnresolvedDrsUrl]] = IO(
    List(
      UnresolvedDrsUrl(fakeDrsUrlWithGcsResolutionOnly, "/path/to/nowhere")
    )
  )

  val fakeAccessInput: IO[List[UnresolvedDrsUrl]] = IO(
    List(
      UnresolvedDrsUrl("https://my-fake-access-url.com", "/path/to/somewhereelse")
    )
  )

  val fakeBulkGoogleInput: IO[List[UnresolvedDrsUrl]] = IO(
    List(
      UnresolvedDrsUrl("drs://my-fake-google-url.com", "/path/to/nowhere"),
      UnresolvedDrsUrl("drs://my-fake-google-url.com2", "/path/to/nowhere2"),
      UnresolvedDrsUrl("drs://my-fake-google-url.com3", "/path/to/nowhere3"),
      UnresolvedDrsUrl("drs://my-fake-google-url.com4", "/path/to/nowhere4")
    )
  )

  val fakeBulkAccessInput: IO[List[UnresolvedDrsUrl]] = IO(
    List(
      UnresolvedDrsUrl("drs://my-fake-access-url.com", "/path/to/somewhereelse"),
      UnresolvedDrsUrl("drs://my-fake-access-url2.com", "/path/to/somewhereelse2"),
      UnresolvedDrsUrl("drs://my-fake-access-url3.com", "/path/to/somewhereelse3"),
      UnresolvedDrsUrl("drs://my-fake-access-url4.com", "/path/to/somewhereelse4")
    )
  )

  behavior of "DrsLocalizerMain"

  it should "fail if drs input is not passed" in {
    DrsLocalizerMain.run(List(fakeDownloadLocation)).unsafeRunSync() shouldBe ExitCode.Error
  }

  it should "fail if download location is not passed" in {
    DrsLocalizerMain.run(List(MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly)).unsafeRunSync() shouldBe ExitCode.Error
  }

  it should "tolerate no URLs being provided" in {
    val mockDownloadFactory = new DownloaderFactory {
      override def buildGcsUriDownloader(gcsPath: String,
                                         serviceAccountJsonOption: Option[String],
                                         downloadLoc: String,
                                         requesterPaysProjectOption: Option[String]
      ): Downloader =
        // This test path should never ask for the Google downloader
        throw new RuntimeException("test failure111")

      override def buildBulkAccessUrlDownloader(urlsToDownload: List[ResolvedDrsUrl]): Downloader =
        // This test path should never ask for the Bulk downloader
        throw new RuntimeException("test failure111")
    }
    val mockdrsLocalizer =
      new MockDrsLocalizerMain(IO(List()), mockDownloadFactory, FakeAccessTokenStrategy, Option(fakeRequesterPaysId))
    val downloaders: List[Downloader] = mockdrsLocalizer.buildDownloaders().unsafeRunSync()
    downloaders.length shouldBe 0
  }

  it should "build correct downloader(s) for a single google URL" in {
    val mockDownloadFactory = new DownloaderFactory {
      override def buildGcsUriDownloader(gcsPath: String,
                                         serviceAccountJsonOption: Option[String],
                                         downloadLoc: String,
                                         requesterPaysProjectOption: Option[String]
      ): Downloader =
        GcsUriDownloader(gcsPath, serviceAccountJsonOption, downloadLoc, requesterPaysProjectOption)

      override def buildBulkAccessUrlDownloader(urlsToDownload: List[ResolvedDrsUrl]): Downloader =
        // This test path should never ask for the Bulk downloader
        throw new RuntimeException("test failure111")
    }

    val mockdrsLocalizer = new MockDrsLocalizerMain(IO(List(fakeGoogleUrls.head._1)),
                                                    mockDownloadFactory,
                                                    FakeAccessTokenStrategy,
                                                    Option(fakeRequesterPaysId)
    )
    val downloaders: List[Downloader] = mockdrsLocalizer.buildDownloaders().unsafeRunSync()
    downloaders.length shouldBe 1

    val correct = downloaders.head match {
      case _: GcsUriDownloader => true
      case _ => false
    }
    correct shouldBe true
  }

  it should "build correct downloader(s) for a single access URL" in {
    val mockDownloadFactory = new DownloaderFactory {
      override def buildGcsUriDownloader(gcsPath: String,
                                         serviceAccountJsonOption: Option[String],
                                         downloadLoc: String,
                                         requesterPaysProjectOption: Option[String]
      ): Downloader =
        // This test path should never ask for the GCS downloader
        throw new RuntimeException("test failure")

      override def buildBulkAccessUrlDownloader(urlsToDownload: List[ResolvedDrsUrl]): Downloader =
        BulkAccessUrlDownloader(urlsToDownload)
    }

    val mockdrsLocalizer = new MockDrsLocalizerMain(IO(List(fakeAccessUrls.head._1)),
                                                    mockDownloadFactory,
                                                    FakeAccessTokenStrategy,
                                                    Option(fakeRequesterPaysId)
    )
    val downloaders: List[Downloader] = mockdrsLocalizer.buildDownloaders().unsafeRunSync()
    downloaders.length shouldBe 1

    val expected = BulkAccessUrlDownloader(
      List(fakeAccessUrls.head._2)
    )
    expected shouldEqual downloaders.head
  }

  it should "build correct downloader(s) for multiple google URLs" in {
    val mockDownloadFactory = new DownloaderFactory {
      override def buildGcsUriDownloader(gcsPath: String,
                                         serviceAccountJsonOption: Option[String],
                                         downloadLoc: String,
                                         requesterPaysProjectOption: Option[String]
      ): Downloader =
        GcsUriDownloader(gcsPath, serviceAccountJsonOption, downloadLoc, requesterPaysProjectOption)

      override def buildBulkAccessUrlDownloader(urlsToDownload: List[ResolvedDrsUrl]): Downloader =
        // This test path should never ask for the GCS downloader
        throw new RuntimeException("test failure")
    }
    val unresolvedUrls: List[UnresolvedDrsUrl] = fakeGoogleUrls.map(pair => pair._1).toList
    val mockdrsLocalizer = new MockDrsLocalizerMain(IO(unresolvedUrls),
                                                    mockDownloadFactory,
                                                    FakeAccessTokenStrategy,
                                                    Option(fakeRequesterPaysId)
    )
    val downloaders: List[Downloader] = mockdrsLocalizer.buildDownloaders().unsafeRunSync()
    downloaders.length shouldBe unresolvedUrls.length

    val countGoogleDownloaders = downloaders.count(downloader =>
      downloader match {
        case _: GcsUriDownloader => true
        case _ => false
      }
    )
    // We expect one GCS downloader for each GCS uri provided
    countGoogleDownloaders shouldBe downloaders.length
  }

  it should "build a single bulk downloader for multiple access URLs" in {
    val mockDownloadFactory = new DownloaderFactory {
      override def buildGcsUriDownloader(gcsPath: String,
                                         serviceAccountJsonOption: Option[String],
                                         downloadLoc: String,
                                         requesterPaysProjectOption: Option[String]
      ): Downloader =
        // This test path should never ask for the GCS downloader
        throw new RuntimeException("test failure")

      override def buildBulkAccessUrlDownloader(urlsToDownload: List[ResolvedDrsUrl]): Downloader =
        BulkAccessUrlDownloader(urlsToDownload)
    }
    val unresolvedUrls: List[UnresolvedDrsUrl] = fakeAccessUrls.map(pair => pair._1).toList
    val mockdrsLocalizer = new MockDrsLocalizerMain(IO(unresolvedUrls),
                                                    mockDownloadFactory,
                                                    FakeAccessTokenStrategy,
                                                    Option(fakeRequesterPaysId)
    )
    val downloaders: List[Downloader] = mockdrsLocalizer.buildDownloaders().unsafeRunSync()
    downloaders.length shouldBe 1

    val countBulkDownloaders = downloaders.count(downloader =>
      downloader match {
        case _: BulkAccessUrlDownloader => true
        case _ => false
      }
    )
    // We expect one total Bulk downloader for all access URIs to share
    countBulkDownloaders shouldBe 1
    val expected = BulkAccessUrlDownloader(
      fakeAccessUrls.map(pair => pair._2).toList
    )
    expected shouldEqual downloaders.head
  }

  it should "build 1 bulk downloader and 5 google downloaders for a mix of URLs" in {
    val unresolvedUrls: List[UnresolvedDrsUrl] =
      fakeAccessUrls.map(pair => pair._1).toList ++ fakeGoogleUrls.map(pair => pair._1).toList
    val mockdrsLocalizer = new MockDrsLocalizerMain(IO(unresolvedUrls),
                                                    DrsLocalizerMain.defaultDownloaderFactory,
                                                    FakeAccessTokenStrategy,
                                                    Option(fakeRequesterPaysId)
    )
    val downloaders: List[Downloader] = mockdrsLocalizer.buildDownloaders().unsafeRunSync()

    downloaders.length shouldBe 6

    // we expect a single bulk downloader despite 5 access URLs being provided
    val countBulkDownloaders = downloaders.count(downloader =>
      downloader match {
        case _: BulkAccessUrlDownloader => true
        case _ => false
      }
    )
    // We expect one GCS downloader for each GCS uri provided
    countBulkDownloaders shouldBe 1
    val countGoogleDownloaders = downloaders.count(downloader =>
      downloader match {
        case _: GcsUriDownloader => true
        case _ => false
      }
    )
    // We expect one GCS downloader for each GCS uri provided
    countBulkDownloaders shouldBe 1
    countGoogleDownloaders shouldBe 5
  }

  it should "accept arguments and run successfully without Requester Pays ID" in {
    val unresolved = fakeGoogleUrls.head._1
    val mockDrsLocalizer = new MockDrsLocalizerMain(IO(List(unresolved)),
                                                    DrsLocalizerMain.defaultDownloaderFactory,
                                                    FakeAccessTokenStrategy,
                                                    None
    )
    val expected = GcsUriDownloader(
      gcsUrl = fakeGoogleUrls.get(unresolved).get.drsResponse.gsUri.get,
      serviceAccountJson = None,
      downloadLoc = unresolved.downloadDestinationPath,
      requesterPaysProjectIdOption = None
    )
    val downloader: Downloader = mockDrsLocalizer.buildDownloaders().unsafeRunSync().head
    downloader shouldBe expected
  }

  it should "run successfully with all 3 arguments" in {
    val unresolved = fakeGoogleUrls.head._1
    val mockDrsLocalizer = new MockDrsLocalizerMain(IO(List(unresolved)),
                                                    DrsLocalizerMain.defaultDownloaderFactory,
                                                    FakeAccessTokenStrategy,
                                                    Option(fakeRequesterPaysId)
    )
    val expected = GcsUriDownloader(
      gcsUrl = fakeGoogleUrls.get(unresolved).get.drsResponse.gsUri.get,
      serviceAccountJson = None,
      downloadLoc = unresolved.downloadDestinationPath,
      requesterPaysProjectIdOption = Option(fakeRequesterPaysId)
    )
    val downloader: Downloader = mockDrsLocalizer.buildDownloaders().unsafeRunSync().head
    downloader shouldBe expected
  }

  it should "successfully identify uri types, preferring access" in {
    val exampleAccessResponse = DrsResolverResponse(accessUrl = Option(AccessUrl("https://something.com", FakeHashes)))
    val exampleGoogleResponse = DrsResolverResponse(gsUri = Option("gs://something"))
    val exampleMixedResponse = DrsResolverResponse(accessUrl = Option(AccessUrl("https://something.com", FakeHashes)),
                                                   gsUri = Option("gs://something")
    )
    DrsLocalizerMain.toValidatedUriType(exampleAccessResponse.accessUrl,
                                        exampleAccessResponse.gsUri
    ) shouldBe URIType.ACCESS
    DrsLocalizerMain.toValidatedUriType(exampleGoogleResponse.accessUrl,
                                        exampleGoogleResponse.gsUri
    ) shouldBe URIType.GCS
    DrsLocalizerMain.toValidatedUriType(exampleMixedResponse.accessUrl,
                                        exampleMixedResponse.gsUri
    ) shouldBe URIType.ACCESS
  }

  it should "throw an exception if the DRS Resolver response is invalid" in {
    val badAccessResponse = DrsResolverResponse(accessUrl = Option(AccessUrl("hQQps://something.com", FakeHashes)))
    val badGoogleResponse = DrsResolverResponse(gsUri = Option("gQQs://something"))
    val emptyResponse = DrsResolverResponse()

    the[RuntimeException] thrownBy {
      DrsLocalizerMain.toValidatedUriType(badAccessResponse.accessUrl, badAccessResponse.gsUri)
    } should have message "Resolved Access URL does not start with https://"

    the[RuntimeException] thrownBy {
      DrsLocalizerMain.toValidatedUriType(badGoogleResponse.accessUrl, badGoogleResponse.gsUri)
    } should have message "Resolved Google URL does not start with gs://"

    the[RuntimeException] thrownBy {
      DrsLocalizerMain.toValidatedUriType(emptyResponse.accessUrl, emptyResponse.gsUri)
    } should have message "DRS response did not contain any URLs"
  }
}

object MockDrsPaths {
  val fakeDrsUrlWithGcsResolutionOnly = "drs://abc/foo-123/abc123"
  val fakeDrsUrlWithAccessUrlResolutionOnly = "drs://def/bar-456/def456"
  val fakeDrsUrlWithAccessUrlAndGcsResolution = "drs://ghi/baz-789/ghi789"
  val fakeDrsUrlWithoutAnyResolution = "drs://foo/bar/no-gcs-path"

  val fakeGoogleUrls: Map[UnresolvedDrsUrl, ResolvedDrsUrl] = Map(
    (UnresolvedDrsUrl("drs://abc/foo-123/google/0", "/path/to/google/local0"),
     ResolvedDrsUrl(DrsResolverResponse(gsUri = Option("gs://some/uri0")), "/path/to/google/local0", URIType.GCS)
    ),
    (UnresolvedDrsUrl("drs://abc/foo-123/google/1", "/path/to/google/local1"),
     ResolvedDrsUrl(DrsResolverResponse(gsUri = Option("gs://some/uri1")), "/path/to/google/local1", URIType.GCS)
    ),
    (UnresolvedDrsUrl("drs://abc/foo-123/google/2", "/path/to/google/local2"),
     ResolvedDrsUrl(DrsResolverResponse(gsUri = Option("gs://some/uri2")), "/path/to/google/local2", URIType.GCS)
    ),
    (UnresolvedDrsUrl("drs://abc/foo-123/google/3", "/path/to/google/local3"),
     ResolvedDrsUrl(DrsResolverResponse(gsUri = Option("gs://some/uri3")), "/path/to/google/local3", URIType.GCS)
    ),
    (UnresolvedDrsUrl("drs://abc/foo-123/google/4", "/path/to/google/local4"),
     ResolvedDrsUrl(DrsResolverResponse(gsUri = Option("gs://some/uri4")), "/path/to/google/local4", URIType.GCS)
    )
  )

  val fakeAccessUrls: Map[UnresolvedDrsUrl, ResolvedDrsUrl] = Map(
    (UnresolvedDrsUrl("drs://abc/foo-123/access/0", "/path/to/access/local0"),
     ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://abc/foo-123/access/0", FakeHashes))),
                    "/path/to/access/local0",
                    URIType.ACCESS
     )
    ),
    (UnresolvedDrsUrl("drs://abc/foo-123/access/1", "/path/to/access/local1"),
     ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://abc/foo-123/access/1", FakeHashes))),
                    "/path/to/access/local1",
                    URIType.ACCESS
     )
    ),
    (UnresolvedDrsUrl("drs://abc/foo-123/access/2", "/path/to/access/local2"),
     ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://abc/foo-123/access/2", FakeHashes))),
                    "/path/to/access/local2",
                    URIType.ACCESS
     )
    ),
    (UnresolvedDrsUrl("drs://abc/foo-123/access/3", "/path/to/access/local3"),
     ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://abc/foo-123/access/3", FakeHashes))),
                    "/path/to/access/local3",
                    URIType.ACCESS
     )
    ),
    (UnresolvedDrsUrl("drs://abc/foo-123/access/4", "/path/to/access/local4"),
     ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://abc/foo-123/access/4", FakeHashes))),
                    "/path/to/access/local4",
                    URIType.ACCESS
     )
    )
  )
}

class MockDrsLocalizerMain(toResolveAndDownload: IO[List[UnresolvedDrsUrl]],
                           downloaderFactory: DownloaderFactory,
                           drsCredentials: DrsCredentials,
                           requesterPaysProjectIdOption: Option[String]
) extends DrsLocalizerMain(toResolveAndDownload,
                           downloaderFactory,
                           FakeAccessTokenStrategy,
                           requesterPaysProjectIdOption
    ) {

  override def getDrsPathResolver: IO[DrsPathResolver] =
    IO {
      new MockDrsLocalizerDrsPathResolver(cloud.nio.impl.drs.MockDrsPaths.mockDrsConfig)
    }
  override def resolveSingleUrl(resolverObject: DrsPathResolver,
                                drsUrlToResolve: UnresolvedDrsUrl
  ): IO[ResolvedDrsUrl] =
    IO {
      if (!fakeAccessUrls.contains(drsUrlToResolve) && !fakeGoogleUrls.contains(drsUrlToResolve)) {
        throw new RuntimeException("Unexpected URI during testing")
      }
      fakeAccessUrls.getOrElse(
        drsUrlToResolve,
        fakeGoogleUrls.getOrElse(drsUrlToResolve, ResolvedDrsUrl(DrsResolverResponse(), "/12/3/", URIType.UNKNOWN))
      )
    }
}

class MockDrsLocalizerDrsPathResolver(drsConfig: DrsConfig)
    extends DrsPathResolver(drsConfig, FakeAccessTokenStrategy) {

  override def resolveDrs(drsPath: String, fields: NonEmptyList[DrsResolverField.Value]): IO[DrsResolverResponse] = {

    val drsResolverResponse = DrsResolverResponse(
      size = Option(1234),
      hashes = FakeHashes
    )

    IO.pure(drsPath) map {
      case MockDrsPaths.fakeDrsUrlWithGcsResolutionOnly =>
        drsResolverResponse.copy(gsUri = Option("gs://abc/foo-123/abc123"))
      case MockDrsPaths.fakeDrsUrlWithoutAnyResolution =>
        drsResolverResponse
      case MockDrsPaths.fakeDrsUrlWithAccessUrlResolutionOnly =>
        drsResolverResponse.copy(accessUrl = Option(AccessUrl(url = "http://abc/def/ghi.bam", headers = None)))
      case MockDrsPaths.fakeDrsUrlWithAccessUrlAndGcsResolution =>
        drsResolverResponse.copy(accessUrl = Option(AccessUrl(url = "http://abc/def/ghi.bam", headers = None)),
                                 gsUri = Option("gs://some/uri")
        )
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
