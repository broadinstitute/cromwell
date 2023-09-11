package drs.localizer.downloaders

import cats.effect.IO
import cloud.nio.impl.drs.{AccessUrl, DrsResolverResponse}
import common.assertion.CromwellTimeoutSpec
import drs.localizer.{ResolvedDrsUrl, URIType}
//import drs.localizer.URIType.URIType
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers


import java.nio.file.Path



class BulkAccessUrlDownloaderSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  val ex1 = ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://my.fake/url123", None))), "path/to/local/download/dest", URIType.HTTPS)
  val ex2 = ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://my.fake/url1234", None))), "path/to/local/download/dest2", URIType.HTTPS)
  val ex3 = ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://my.fake/url1235", None))), "path/to/local/download/dest3", URIType.HTTPS)
  val emptyList : List[ResolvedDrsUrl] = List()
  val oneElement: List[ResolvedDrsUrl] = List(ex1)
  val threeElements: List[ResolvedDrsUrl] = List(ex1, ex2, ex3)

  it should "correctly parse a collection of Access Urls into a manifest.json" in {
    val expected: String =
      s"""|[
          |  {
          |    "url" : "https://my.fake/url123",
          |    "filepath" : "path/to/local/download/dest"
          |  },
          |  {
          |    "url" : "https://my.fake/url1234",
          |    "filepath" : "path/to/local/download/dest2"
          |  },
          |  {
          |    "url" : "https://my.fake/url1235",
          |    "filepath" : "path/to/local/download/dest3"
          |  }
          |]""".stripMargin

    val downloader = BulkAccessUrlDownloader(threeElements)

    val filepath: IO[Path] = downloader.generateJsonManifest(threeElements)
    val source = scala.io.Source.fromFile(filepath.unsafeRunSync().toString)
    val lines = try source.mkString finally source.close()
    lines shouldBe expected
  }

  it should "properly construct empty JSON array from empty list." in {
    val expected: String =
      s"""|[
          |
          |]""".stripMargin

    val downloader = BulkAccessUrlDownloader(emptyList)
    val filepath: IO[Path] = downloader.generateJsonManifest(emptyList)
    val source = scala.io.Source.fromFile(filepath.unsafeRunSync().toString)
    val lines = try source.mkString finally source.close()
    lines shouldBe expected
  }

  it should "properly construct JSON array from single element list." in {
    val expected: String =
      s"""|[
          |  {
          |    "url" : "https://my.fake/url123",
          |    "filepath" : "path/to/local/download/dest"
          |  }
          |]""".stripMargin

    val downloader = BulkAccessUrlDownloader(oneElement)
    val filepath: IO[Path] = downloader.generateJsonManifest(oneElement)
    val source = scala.io.Source.fromFile(filepath.unsafeRunSync().toString)
    val lines = try source.mkString finally source.close()
    lines shouldBe expected
  }

  it should "properly construct the invocation command" in {
    val downloader = BulkAccessUrlDownloader(oneElement)
    val filepath: Path = downloader.generateJsonManifest(threeElements).unsafeRunSync()
    val expected = s"""getm --manifest ${filepath.toString}"""
    downloader.generateGetmCommand(filepath) shouldBe expected
  }
}
