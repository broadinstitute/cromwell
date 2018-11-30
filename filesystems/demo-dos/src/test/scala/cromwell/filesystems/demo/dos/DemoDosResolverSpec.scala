//package cromwell.filesystems.demo.dos
//
//import java.util.UUID
//
//import com.typesafe.config.{Config, ConfigFactory}
//import cromwell.filesystems.demo.dos.DemoDosResolverSpec._
//import org.apache.http._
//import org.apache.http.client.methods.{CloseableHttpResponse, HttpPost}
//import org.apache.http.entity.{ContentType, StringEntity}
//import org.apache.http.impl.EnglishReasonPhraseCatalog
//import org.apache.http.impl.client.CloseableHttpClient
//import org.apache.http.message.{BasicHttpResponse, BasicStatusLine}
//import org.apache.http.protocol.HttpContext
//import org.apache.http.util.EntityUtils
//import org.scalatest.prop.TableDrivenPropertyChecks
//import org.scalatest.{FlatSpec, Matchers}
//
//class DemoDosResolverSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {
//
//  behavior of "DemoDosResolver"
//
//  private val demoDosPathBuilder = new DemoDosPathBuilder()
//
//  private val jsonTemplate = s"""{\\"dosUrl\\": \\"$${dosPath}\\"}"""
//
//  private val exampleConfig = ConfigFactory.parseString(
//    s"""|demo.dos.martha.url = "http://example.org/some/possible/path"
//        |demo.dos.martha.request.json-template = "$jsonTemplate"
//        |demo.dos.martha.response.jq-filter = ".dos"
//        |""".stripMargin
//  )
//
//  private val debugConfig = ConfigFactory.parseString(
//    s"""|demo.dos.martha.url = "http://example.org/some/possible/path"
//        |demo.dos.martha.request.json-template = "$jsonTemplate"
//        |demo.dos.martha.response.jq-filter = ".dos"
//        |demo.dos.martha.debug = true
//        |""".stripMargin
//  )
//
//  private val jqGoogleFilter = """.dos.data_object.urls[] | select(.url | startswith(\"gs://\")) | .url"""
//  private val filterConfig = ConfigFactory.parseString(
//    s"""|demo.dos.martha.url = "http://example.org/some/possible/path"
//        |demo.dos.martha.request.json-template = "$jsonTemplate"
//        |demo.dos.martha.response.jq-filter = "$jqGoogleFilter"
//        |""".stripMargin
//  )
//
//  private val okResponseStatusLine = statusLine(HttpStatus.SC_OK)
//
//  private val forbiddenResponseStatusLine = statusLine(HttpStatus.SC_FORBIDDEN)
//
//  private val stringResponseEntity = new StringEntity(
//    """|{
//       |  "dos": "gs://my-gs-bucket/path/to/file.txt"
//       |}
//       |""".stripMargin,
//    ContentType.APPLICATION_JSON
//  )
//
//  private val arrayResponseEntity = new StringEntity(
//    """|{
//       |  "dos": [
//       |    "s3://my-s3-bucket/path/to/file.txt",
//       |    "gs://my-gs-bucket/path/to/file.txt",
//       |    "oss://my-oss-bucket/path/to/file.txt"
//       |  ]
//       |}
//       |""".stripMargin,
//    ContentType.APPLICATION_JSON
//  )
//
//  private val gcsFirstResponseEntity = new StringEntity(
//    """|{
//       |  "dos": {
//       |    "data_object": {
//       |      "checksums": [
//       |        {
//       |          "checksum": "559573fbc52b3b3ca71eaea7fb22be5009a6bc54",
//       |          "type": "md5"
//       |        }
//       |      ],
//       |      "description": "",
//       |      "id": "630d31c3-381e-488d-b639-ce5d047a0142",
//       |      "mime_type": "",
//       |      "size": 2201638,
//       |      "urls": [
//       |        {
//       |          "url": "gs://my-gs-bucket/path/to/file.txt"
//       |        },
//       |        {
//       |          "url": "s3://my-s3-bucket/path/to/file.txt"
//       |        }
//       |      ],
//       |      "version": "2018-05-26T134315.070672Z"
//       |    }
//       |  },
//       |  "googleServiceAccount": {
//       |  }
//       |}
//       |""".stripMargin,
//    ContentType.APPLICATION_JSON
//  )
//
//  private val gcsSecondResponseEntity = new StringEntity(
//    """|{
//       |  "dos": {
//       |    "data_object": {
//       |      "checksums": [
//       |        {
//       |          "checksum": "fdd7e5dc38915cc0a1f72a97cb9195189b4098eb",
//       |          "type": "md5"
//       |        }
//       |      ],
//       |      "description": "",
//       |      "id": "01b048d0-e128-4cb0-94e9-b2d2cab7563d",
//       |      "mime_type": "",
//       |      "size": 37501686827,
//       |      "urls": [
//       |        {
//       |          "url": "s3://my-s3-bucket/path/to/file.txt"
//       |        },
//       |        {
//       |          "url": "gs://my-gs-bucket/path/to/file.txt"
//       |        }
//       |      ],
//       |      "version": "2018-05-26T133719.491781Z"
//       |    }
//       |  },
//       |  "googleServiceAccount": {
//       |  }
//       |}
//       |""".stripMargin,
//    ContentType.APPLICATION_JSON
//  )
//
//  private val parseTests = Table(
//    ("description", "config", "responseEntity"),
//    ("find gs paths from string json responses", exampleConfig, stringResponseEntity),
//    ("find gs paths from array json responses", exampleConfig, arrayResponseEntity),
//    ("find gs paths when first in an array", filterConfig, gcsFirstResponseEntity),
//    ("find gs paths when second in an array", filterConfig, gcsSecondResponseEntity),
//  )
//
//  forAll(parseTests) { (description, config, responseEntity) =>
//    it should description in {
//      val mockClient = new MockClient(okResponseStatusLine, Option(responseEntity))
//
//      val resolver = new MockDemoDosResolver(config, mockClient)
//
//      val uuid = UUID.randomUUID()
//      val dosUrl = s"dos://my-dos-host/$uuid"
//      val demoDosPath = demoDosPathBuilder.build(dosUrl).get.asInstanceOf[DemoDosPath]
//      resolver.getContainerRelativePath(demoDosPath) should be("my-gs-bucket/path/to/file.txt")
//      mockClient.targetOption.get should be(new HttpHost("example.org", -1, "http"))
//      mockClient.contextOption should be(empty)
//
//      mockClient.requestOption.get should be(an[HttpPost])
//      val postRequest = mockClient.requestOption.get.asInstanceOf[HttpPost]
//      postRequest.getEntity.getContentType.getValue should be(ContentType.APPLICATION_JSON.toString)
//      EntityUtils.toString(postRequest.getEntity) should be(s"""{"dosUrl": "$dosUrl"}""")
//    }
//  }
//
//  private val normalExceptionMessage =
//    s"""Unexpected response looking up $${dosUrl} from http://example.org/some/possible/path."""
//
//  private val debugExceptionMessage =
//    s"""|Unexpected response looking up $${dosUrl} from http://example.org/some/possible/path.
//        |HTTP/1.1 403 Forbidden []""".stripMargin
//
//
//  private val forbiddenTests = Table(
//    ("description", "config", "exceptionMessageTemplate"),
//    ("throw an exception for a forbidden request", exampleConfig, normalExceptionMessage),
//    ("throw an exception for a debugged forbidden request", debugConfig, debugExceptionMessage)
//  )
//
//  forAll(forbiddenTests) { (description, config, exceptionMessage) =>
//    it should description in {
//      val mockClient = new MockClient(forbiddenResponseStatusLine)
//
//      val resolver = new MockDemoDosResolver(config, mockClient)
//
//      val uuid = UUID.randomUUID()
//      val dosUrl = s"dos://my-dos-host/$uuid"
//      val demoDosPath = demoDosPathBuilder.build(dosUrl).get.asInstanceOf[DemoDosPath]
//      val expectedMessage = exceptionMessage.replace(s"$${dosUrl}", dosUrl)
//      the[RuntimeException] thrownBy resolver.getContainerRelativePath(demoDosPath) should have message expectedMessage
//    }
//  }
//}
//
//object DemoDosResolverSpec {
//
//  private def statusLine(status: Int) = new BasicStatusLine(
//    HttpVersion.HTTP_1_1,
//    status,
//    EnglishReasonPhraseCatalog.INSTANCE.getReason(status, null)
//  )
//
//  class MockDemoDosResolver(config: Config, mockClient: MockClient)
//    extends DemoDosResolver(config) {
//      override def createClient(): CloseableHttpClient = mockClient
//  }
//
//  class MockClient(responseStatusLine: StatusLine, responseEntityOption: Option[HttpEntity] = None)
//    extends CloseableHttpClient {
//    var targetOption: Option[HttpHost] = None
//    var requestOption: Option[HttpRequest] = None
//    var contextOption: Option[HttpContext] = None
//    var responseOption: Option[HttpResponse] = None
//
//    override def doExecute(target: HttpHost, request: HttpRequest, context: HttpContext): CloseableHttpResponse = {
//      targetOption = Option(target)
//      requestOption = Option(request)
//      contextOption = Option(context)
//      val response = new MockResponse(responseStatusLine)
//      responseEntityOption.foreach(response.setEntity)
//      responseOption = Option(response)
//      response
//    }
//
//    override def close(): Unit = {}
//
//    override def getParams = throw new UnsupportedOperationException()
//
//    override def getConnectionManager = throw new UnsupportedOperationException()
//  }
//
//  class MockResponse(statusLine: StatusLine) extends BasicHttpResponse(statusLine) with CloseableHttpResponse {
//    override def close(): Unit = {}
//  }
//
//}
