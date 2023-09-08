package drs.localizer.downloaders

import cats.effect.IO
import cloud.nio.impl.drs.DrsResolverField.AccessUrl
import cloud.nio.impl.drs.{AccessUrl, DrsResolverResponse}
import common.assertion.CromwellTimeoutSpec
import drs.localizer.{ResolvedDrsUrl, URIType}
//import drs.localizer.URIType.URIType
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.Path



class BulkAccessUrlDownloaderSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {
  val emptyList : List[ResolvedDrsUrl] = List()
  val oneElement: List[ResolvedDrsUrl] = List(ResolvedDrsUrl(DrsResolverResponse(), "path/to/local/download/dest", URIType.HTTPS))
  val threeElements: List[ResolvedDrsUrl] = List(
    ResolvedDrsUrl(DrsResolverResponse(accessUrl = Option(AccessUrl("https://my.fake/url123", None))), "path/to/local/download/dest", URIType.HTTPS),
    ResolvedDrsUrl(DrsResolverResponse(), "path/to/local/download/dest2", URIType.HTTPS),
    ResolvedDrsUrl(DrsResolverResponse(), "path/to/local/download/dest3", URIType.HTTPS))

  it should "correctly parse a collection of Access Urls into a manifest.json" in {
    val expected : String =
      s"""[
         |  {
         |    "url" : "https://my.fake/url123",
         |    "filepath" : "./test1.out",
         |    "checksum" : "Valid()",
         |    "checksum-algorithm" : "null"
         |  },
         |  {
         |    "url" : "https://my.fake/url1234",
         |    "filepath" : "./test2.out",
         |    "checksum" : "Valid()",
         |    "checksum-algorithm" : "null"
         |  }
         |]""".stripMargin

    val downloader = BulkAccessUrlDownloader(threeElements)

    val filepath : IO[Path] = downloader.generateJsonManifest(threeElements)
    val source = scala.io.Source.fromFile(filepath.unsafeRunSync().toString)
    val lines = try source.mkString finally source.close()

    lines shouldBe expected
  }

}
