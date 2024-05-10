package drs.localizer.downloaders

import cats.effect.{ExitCode, IO}
import cloud.nio.impl.drs.{AccessUrl, DrsResolverResponse}
import common.assertion.CromwellTimeoutSpec
import org.scalatest.prop.TableDrivenPropertyChecks._
import drs.localizer.{ResolvedDrsUrl, URIType}

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.Path

class BulkAccessUrlDownloaderSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  val ex1 = ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://my.fake/url123", None))),
                           "path/to/local/download/dest",
                           URIType.ACCESS
  )
  val ex2 = ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://my.fake/url1234", None))),
                           "path/to/local/download/dest2",
                           URIType.ACCESS
  )
  val ex3 = ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://my.fake/url1235", None))),
                           "path/to/local/download/dest3",
                           URIType.ACCESS
  )
  val emptyList: List[ResolvedDrsUrl] = List()
  val oneElement: List[ResolvedDrsUrl] = List(ex1)
  val threeElements: List[ResolvedDrsUrl] = List(ex1, ex2, ex3)

  it should "correctly parse a collection of Access Urls into a manifest.json" in {
    val expected =
      """[{
        |  "url": "https://my.fake/url123",
        |  "filepath": "path/to/local/download/dest"
        |}, {
        |  "url": "https://my.fake/url1234",
        |  "filepath": "path/to/local/download/dest2"
        |}, {
        |  "url": "https://my.fake/url1235",
        |  "filepath": "path/to/local/download/dest3"
        |}]""".stripMargin
    val downloader = BulkAccessUrlDownloader(threeElements)

    val filepath: IO[Path] = downloader.generateJsonManifest(threeElements)
    val source = scala.io.Source.fromFile(filepath.unsafeRunSync().toString)
    val lines =
      try source.mkString
      finally source.close()
    lines shouldBe expected
  }

  it should "properly construct empty JSON array from empty list." in {
    val expected: String = "[]"
    val downloader = BulkAccessUrlDownloader(emptyList)
    val filepath: IO[Path] = downloader.generateJsonManifest(emptyList)
    val source = scala.io.Source.fromFile(filepath.unsafeRunSync().toString)
    val lines =
      try source.mkString
      finally source.close()
    lines shouldBe expected
  }

  it should "properly construct JSON array from single element list." in {
    val expected: String =
      s"""|[{
          |  "url": "https://my.fake/url123",
          |  "filepath": "path/to/local/download/dest"
          |}]""".stripMargin

    val downloader = BulkAccessUrlDownloader(oneElement)
    val filepath: IO[Path] = downloader.generateJsonManifest(oneElement)
    val source = scala.io.Source.fromFile(filepath.unsafeRunSync().toString)
    val lines =
      try source.mkString
      finally source.close()
    lines shouldBe expected
  }

  it should "properly construct the invocation command" in {
    val downloader = BulkAccessUrlDownloader(oneElement)
    val filepath: Path = downloader.generateJsonManifest(threeElements).unsafeRunSync()
    val expected = s"""timeout 24h getm --manifest ${filepath.toString} -vv"""
    downloader.generateGetmCommand(filepath) shouldBe expected
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
      (0,
       """ERROR:getm.cli possibly some words "status_code": 503 words""",
       UnrecognizedRetryableDownloadFailure(ExitCode(0))
      ),
      // Recognized because of non-zero exit status and an HTTP status.
      (1,
       """ERROR:getm.cli possibly some words "status_code": 503 words""",
       RecognizedRetryableDownloadFailure(ExitCode(1))
      ),
      // Recognized because of non-zero exit status and an HTTP status.
      (1,
       """ERROR:getm.cli possibly some words "status_code": 408 more words""",
       RecognizedRetryableDownloadFailure(ExitCode(1))
      ),
      // Recognized and non-retryable because of non-zero exit status and 404 HTTP status.
      (1,
       """ERROR:getm.cli possibly some words "status_code": 404 even more words""",
       FatalDownloadFailure(ExitCode(1))
      ),
      // Unrecognized because of zero exit status and 404 HTTP status.
      (0,
       """ERROR:getm.cli possibly some words "status_code": 404 even more words""",
       UnrecognizedRetryableDownloadFailure(ExitCode(0))
      )
    )
    val bulkDownloader = BulkAccessUrlDownloader(null)

    forAll(results) { (exitCode, stderr, expected) =>
      it should s"produce $expected for exitCode $exitCode and stderr '$stderr'" in {
        bulkDownloader.toDownloadResult(GetmResult(exitCode, stderr)) shouldBe expected
      }
    }
  }
}
