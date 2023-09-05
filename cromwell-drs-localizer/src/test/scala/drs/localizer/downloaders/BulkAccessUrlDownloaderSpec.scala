package drs.localizer.downloaders

import cloud.nio.impl.drs.AccessUrl
import common.assertion.CromwellTimeoutSpec
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.Path
import scala.util.Try


class BulkAccessUrlDownloaderSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers {

  it should "correctly parse a collection of Access Urls into a manifest.json" in {

    val fakeAccessURLToDownloadLocationMap : Map[AccessUrl, String] = Map(
      (AccessUrl("https://my.fake/url123", None) -> "./test1.out"),
      (AccessUrl("https://my.fake/url1234", None) -> "./test2.out")
    )

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

    val downloader = BulkAccessUrlDownloader(fakeAccessURLToDownloadLocationMap)

    val filepath : Try[Path] = downloader.generateJsonManifest(fakeAccessURLToDownloadLocationMap)
    val source = scala.io.Source.fromFile(filepath.get.toString)
    val lines = try source.mkString finally source.close()

    lines shouldBe expected
  }
}
