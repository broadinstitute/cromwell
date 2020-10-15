package cloud.nio.impl.drs

import java.nio.channels.ReadableByteChannel
import java.time.{Instant, OffsetDateTime, ZoneOffset}

import cats.data.NonEmptyList
import cats.effect.IO
import cloud.nio.impl.drs.DrsCloudNioFileProvider.DrsReadInterpreter
import cloud.nio.spi.CloudNioFileList
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
    fileSystemProvider.drsConfig should be(DrsConfig("https://from.config"))
    fileSystemProvider.accessTokenAcceptableTTL should be(1.minute)
    fileSystemProvider.fileProvider should be(a[DrsCloudNioFileProvider])
    fileSystemProvider.isFatal(new RuntimeException) should be(false)
    fileSystemProvider.isTransient(new RuntimeException) should be(false)
    fileSystemProvider.getScheme should be("drs")
    fileSystemProvider.getHost("drs://dg.123/abc") should be("dg.123")
    fileSystemProvider.getHost("drs://dg.example.com/abc") should be("dg.example.com")
  }

  it should "be able to get the hostname from variously formatted DRS URIs" in {
    val config = ConfigFactory.parseString(
      """martha.url = "https://from.config"
        |access-token-acceptable-ttl = 1 minute
        |""".stripMargin
    )

    val fileSystemProvider = new MockDrsCloudNioFileSystemProvider(config = config)
    fileSystemProvider.drsConfig should be(DrsConfig("https://from.config"))
    fileSystemProvider.accessTokenAcceptableTTL should be(1.minute)
    fileSystemProvider.fileProvider should be(a[DrsCloudNioFileProvider])
    fileSystemProvider.isFatal(new RuntimeException) should be(false)
    fileSystemProvider.isTransient(new RuntimeException) should be(false)
    fileSystemProvider.getScheme should be("drs")
    fileSystemProvider.getHost("drs://dg.4503:dg.4503/abc") should be("dg.4503")
    fileSystemProvider.getHost("drs://dg.4503:abc") should be("dg.4503")
    fileSystemProvider.getHost("drs://dg.4503:") should be("dg.4503")
    fileSystemProvider.getHost("drs://dg.4DFC:abc") should be("dg.4DFC")
    fileSystemProvider.getHost("drs://dg.712c:abc") should be("dg.712c")
    fileSystemProvider.getHost("drs://dg.ANV0:abc") should be("dg.ANV0")
    fileSystemProvider.getHost("drs://dg.F82A1A:abc") should be("dg.F82A1A")
    fileSystemProvider.getHost("drs://:abc") should be("")
  }

  it should "list existing drs objects" in {
    val httpResponse = mock[CloseableHttpResponse].smart
    httpResponse.getStatusLine returns new BasicStatusLine(HttpVersion.HTTP_1_1, 200, "OK")

    val httpClient = spy(new MockHttpClient()).smart
    httpClient.setResponseCode(HttpStatus.SC_OK)
    doReturn(httpResponse).when(httpClient).execute(anyObject[HttpPost])

    val httpClientBuilder = mock[HttpClientBuilder].smart
    httpClientBuilder.build() returns httpClient

    val fileSystemProvider = new MockDrsCloudNioFileSystemProvider(httpClientBuilder = httpClientBuilder)
    val fileProvider = fileSystemProvider.fileProvider.asInstanceOf[DrsCloudNioFileProvider]
    fileProvider.existsPaths("dg.123", "abc") should be(true)
    fileProvider.listObjects("dg.123", "abc", None) should be(CloudNioFileList(List("abc"), None))
    fileProvider.listObjects("dg.123", "abc/", None) should be(CloudNioFileList(Nil, None))
  }

  it should "return a file provider that can read bytes" in {
    val drsPathResolver = new MockEngineDrsPathResolver() {
      override def resolveDrsThroughMartha(drsPath: String,
                                           fields: NonEmptyList[MarthaField.Value],
                                          ): IO[MarthaResponse] = {
        IO(
          MarthaResponse(
            size = None,
            timeCreated = None,
            timeUpdated = None,
            gsUri = Option("gs://bucket/object/path"),
            googleServiceAccount = None,
            fileName = None,
            hashes = None,
          )
        )
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
            gsUri = None,
            googleServiceAccount = None,
            fileName = None,
            hashes = Option(Map("rot13" -> "gg0217869")),
          )
        )
      }
    }

    val fileSystemProvider = new MockDrsCloudNioFileSystemProvider(mockResolver = Option(drsPathResolver))
    val fileAttributes = fileSystemProvider.fileProvider.fileAttributes("dg.123", "abc")
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
  }
}
