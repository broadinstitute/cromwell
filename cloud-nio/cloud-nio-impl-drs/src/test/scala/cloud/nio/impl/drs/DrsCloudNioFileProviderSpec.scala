package cloud.nio.impl.drs

import java.nio.channels.ReadableByteChannel
import java.time.{Instant, OffsetDateTime, ZoneOffset}

import cats.data.NonEmptyList
import cats.effect.IO
import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import com.google.api.client.testing.http.apache.MockHttpClient
import com.typesafe.config.ConfigFactory
import common.assertion.CromwellTimeoutSpec
import org.apache.commons.compress.utils.SeekableInMemoryByteChannel
import org.apache.http.client.methods.{CloseableHttpResponse, HttpPost}
import org.apache.http.impl.client.HttpClientBuilder
import org.apache.http.message.BasicStatusLine
import org.apache.http.{HttpStatus, HttpVersion}
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito

import scala.concurrent.duration._

class DrsCloudNioFileProviderSpec extends AnyFlatSpecLike with CromwellTimeoutSpec with Matchers with Mockito {
  behavior of "DrsCloudNioFileProvider"

  it should "parse a config and create a working file system provider" in {
    val config = ConfigFactory.parseString(
      """martha.url = "https://from.config"
        |access-token-acceptable-ttl = 1 minute
        |""".stripMargin
    )

    val fileSystemProvider = new MockDrsCloudNioFileSystemProvider(config = config)
    fileSystemProvider.drsConfig.marthaUrl should be("https://from.config")
    fileSystemProvider.accessTokenAcceptableTTL should be(1.minute)
    fileSystemProvider.fileProvider should be(a[DrsCloudNioFileProvider])
    fileSystemProvider.isFatal(new RuntimeException) should be(false)
    fileSystemProvider.isTransient(new RuntimeException) should be(false)
    fileSystemProvider.getScheme should be("drs")
    val exampleUris = List(
      "drs://dg.123/abc",
      "drs://dg.example.com/abc",
      "drs://drs.data.humancellatlas.org/8aca942c-17f7-4e34-b8fd-3c12e50f9291?version=2019-07-04T151444.185805Z",
      "drs://jade.datarepo-dev.broadinstitute.org" +
        "/v1_93dc1e76-8f1c-4949-8f9b-07a087f3ce7b_8b07563a-542f-4b5c-9e00-e8fe6b1861de",
      "drs://dg.4503/fc046e84-6cf9-43a3-99cc-ffa2964b88cb",
      // Example non-W3C/IETF URIs provided by
      // https://docs.google.com/document/d/1Wf4enSGOEXD5_AE-uzLoYqjIp5MnePbZ6kYTVFp1WoM/edit
      "drs://dg.4DFC:0027045b-9ed6-45af-a68e-f55037b5184c",
      "drs://dg.4503:dg.4503/fc046e84-6cf9-43a3-99cc-ffa2964b88cb",
      "drs://dg.ANV0:dg.ANV0/0db6577e-57bd-48a1-93c6-327c292bcb6b",
      "drs://dg.F82A1A:ed6be7ab-068e-46c8-824a-f39cfbb885cc",
    )
    for (exampleUri <- exampleUris) {
      fileSystemProvider.getHost(exampleUri) should be(exampleUri)
    }
  }

  it should "check existing drs objects" in {
    val httpResponse = mock[CloseableHttpResponse].smart
    httpResponse.getStatusLine returns new BasicStatusLine(HttpVersion.HTTP_1_1, 200, "OK")

    val httpClient = spy(new MockHttpClient()).smart
    httpClient.setResponseCode(HttpStatus.SC_OK)
    doReturn(httpResponse).when(httpClient).execute(anyObject[HttpPost])

    val httpClientBuilder = mock[HttpClientBuilder].smart
    httpClientBuilder.build() returns httpClient

    val fileSystemProvider = new MockDrsCloudNioFileSystemProvider(httpClientBuilder = Option(httpClientBuilder))
    val fileProvider = fileSystemProvider.fileProvider.asInstanceOf[DrsCloudNioFileProvider]
    fileProvider.existsPaths("drs://dg.123/abc", "") should be(true)
  }

  it should "return a file provider that can read bytes" in {
    val drsPathResolver = new MockEngineDrsPathResolver() {
      override def resolveDrsThroughMartha(drsPath: String,
                                           fields: NonEmptyList[MarthaField.Value],
                                          ): IO[MarthaResponse] = {
        IO(MarthaResponse(gsUri = Option("gs://bucket/object/path")))
      }
    }

    val readChannel: ReadableByteChannel = new SeekableInMemoryByteChannel(Array.emptyByteArray)
    val drsReadInterpreter: DrsReadInterpreter = (gsUriOption, googleServiceAccountOption) => {
      IO(
        (gsUriOption, googleServiceAccountOption) match {
          case (Some("gs://bucket/object/path"), None) => readChannel
          case _ => fail(s"Unexpected parameters passed: ($gsUriOption,$googleServiceAccountOption")
        }
      )
    }

    val fileSystemProvider = new MockDrsCloudNioFileSystemProvider(
      mockResolver = Option(drsPathResolver),
      drsReadInterpreter = drsReadInterpreter,
    )
    fileSystemProvider.fileProvider.read("dg.123", "abc", 0) should be(readChannel)
  }

  it should "return a file provider that can return file attributes" in {
    val drsPathResolver = new MockEngineDrsPathResolver() {
      override def resolveDrsThroughMartha(drsPath: String,
                                           fields: NonEmptyList[MarthaField.Value],
                                          ): IO[MarthaResponse] = {
        val instantCreated = Instant.ofEpochMilli(123L)
        val instantUpdated = Instant.ofEpochMilli(456L)
        IO(
          MarthaResponse(
            size = Option(789L),
            timeCreated = Option(OffsetDateTime.ofInstant(instantCreated, ZoneOffset.UTC).toString),
            timeUpdated = Option(OffsetDateTime.ofInstant(instantUpdated, ZoneOffset.UTC).toString),
            hashes = Option(Map("rot13" -> "gg0217869")),
          )
        )
      }
    }

    val fileSystemProvider = new MockDrsCloudNioFileSystemProvider(mockResolver = Option(drsPathResolver))
    val fileAttributes = fileSystemProvider.fileProvider.fileAttributes("drs://dg.123/abc", "")
    val drsFileAttributes = fileAttributes.get.asInstanceOf[DrsCloudNioRegularFileAttributes]
    drsFileAttributes.fileKey() should be("drs://dg.123/abc")
    drsFileAttributes.creationTime().toMillis should be(123L)
    drsFileAttributes.lastModifiedTime().toMillis should be(456L)
    drsFileAttributes.size() should be(789L)
    drsFileAttributes.fileHash should be(Option("gg0217869"))
  }

  it should "throw exceptions for unsupported methods" in {
    val fileSystemProvider = new MockDrsCloudNioFileSystemProvider()
    val fileProvider = fileSystemProvider.fileProvider.asInstanceOf[DrsCloudNioFileProvider]

    the[UnsupportedOperationException] thrownBy {
      fileProvider.copy("", "", "", "")
    } should have message "DRS currently doesn't support copy."

    the[UnsupportedOperationException] thrownBy {
      fileProvider.write("", "")
    } should have message "DRS currently doesn't support write."

    the[UnsupportedOperationException] thrownBy {
      fileProvider.deleteIfExists("", "")
    } should have message "DRS currently doesn't support delete."

    the[UnsupportedOperationException] thrownBy {
      fileProvider.listObjects("", "", None)
    } should have message "DRS currently doesn't support list."
  }
}
