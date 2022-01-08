package cromwell.filesystems.drs

import cloud.nio.impl.drs.{AccessUrl, DrsPathResolver, MarthaResponse, MockEngineDrsPathResolver}
import common.assertion.CromwellTimeoutSpec
import cromwell.cloudsupport.gcp.auth.MockAuthMode
import cromwell.core.WorkflowOptions
import org.apache.http.client.methods.{CloseableHttpResponse, HttpGet}
import org.apache.http.entity.ByteArrayEntity
import org.apache.http.impl.client.{CloseableHttpClient, HttpClientBuilder}
import org.mockito.Mockito.verify
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito
import org.testcontainers.shaded.org.apache.commons.io.IOUtils

import java.nio.ByteBuffer

class DrsReaderSpec extends AnyFlatSpecLike with CromwellTimeoutSpec with Matchers with Mockito {

  behavior of "DrsReader"

  it should "return a gcs reader for a gcs url" in {
    val googleAuthMode = MockAuthMode("unused")
    val workflowOptions = WorkflowOptions.empty
    val requesterPaysProjectIdOption = None
    val drsPathResolver = new MockEngineDrsPathResolver()
    val gsUri = "gs://bucket/object"
    val googleServiceAccount = None
    val marthaResponse = MarthaResponse(gsUri = Option(gsUri), googleServiceAccount = googleServiceAccount)
    val readerIo =
      DrsReader.reader(Option(googleAuthMode), workflowOptions, requesterPaysProjectIdOption, drsPathResolver, marthaResponse)
    readerIo.unsafeRunSync() should be(
      GcsReader(googleAuthMode, workflowOptions, requesterPaysProjectIdOption, gsUri, googleServiceAccount)
    )
  }

  it should "return an access url reader for an access url" in {
    val googleAuthMode = MockAuthMode("unused")
    val workflowOptions = WorkflowOptions.empty
    val requesterPaysProjectIdOption = None
    val drsPathResolver = new MockEngineDrsPathResolver()
    val accessUrl = AccessUrl("https://host/object/path", Option(Map("hello" -> "world")))
    val marthaResponse = MarthaResponse(accessUrl = Option(accessUrl))
    val readerIo =
      DrsReader.reader(Option(googleAuthMode), workflowOptions, requesterPaysProjectIdOption, drsPathResolver, marthaResponse)
    readerIo.unsafeRunSync() should be(
      AccessUrlReader(drsPathResolver, accessUrl)
    )
  }

  it should "return an error without a gcs nor an access url" in {
    val googleAuthMode = MockAuthMode("unused")
    val workflowOptions = WorkflowOptions.empty
    val requesterPaysProjectIdOption = None
    val drsPathResolver = new MockEngineDrsPathResolver()
    val marthaResponse = MarthaResponse()
    val readerIo =
      DrsReader.reader(Option(googleAuthMode), workflowOptions, requesterPaysProjectIdOption, drsPathResolver, marthaResponse)
    the[RuntimeException] thrownBy {
      readerIo.unsafeRunSync()
    } should have message DrsPathResolver.ExtractUriErrorMsg
  }

  it should "return a closeable channel for an access url" in {
    val exampleBytes = Array[Byte](1, 2, 3)
    val httpResponse = mock[CloseableHttpResponse].smart
    httpResponse.getEntity returns new ByteArrayEntity(exampleBytes)

    val httpClient = mock[CloseableHttpClient].smart
    doReturn(httpResponse).when(httpClient).execute(anyObject[HttpGet])

    val httpClientBuilder = mock[HttpClientBuilder].smart
    httpClientBuilder.build() returns httpClient

    val drsPathResolver = new MockEngineDrsPathResolver(httpClientBuilderOverride = Option(httpClientBuilder))

    val accessUrl = AccessUrl("https://host/object/path", Option(Map("hello" -> "world")))
    val marthaResponse = MarthaResponse(accessUrl = Option(accessUrl))
    val channelIo =
      DrsReader.readInterpreter(Option(MockAuthMode("unused")), WorkflowOptions.empty, None)(drsPathResolver, marthaResponse)
    val channel = channelIo.unsafeRunSync()

    val buffer = ByteBuffer.allocate(exampleBytes.length)
    channel.isOpen should be (true)
    IOUtils.readFully(channel, buffer)
    channel.close()

    val httpGetCapture = capture[HttpGet]
    channel.isOpen should be (false)
    buffer.array() should be(exampleBytes)
    verify(httpClient).execute(httpGetCapture.capture)
    verify(httpClient).close()
    verify(httpResponse).close()

    val actualHeaders = httpGetCapture.value.getAllHeaders
    actualHeaders.length should be(1)
    actualHeaders(0).getName should be("hello")
    actualHeaders(0).getValue should be("world")
  }
}
