package cromwell.filesystems.drs

import java.nio.channels.ReadableByteChannel

import cats.effect.IO
import cloud.nio.impl.drs.MarthaResponse
import com.google.cloud.NoCredentials
import com.typesafe.config.{Config, ConfigFactory}
import org.apache.http.impl.client.HttpClientBuilder
import org.scalatest.{FlatSpec, Matchers}

import scala.collection.JavaConverters._


class DrsResolverSpec extends FlatSpec with Matchers {

  private val marthaV2Config: Config = ConfigFactory.parseMap(
    Map(
      "martha.url" -> "https://martha-url/martha_v2",
      "martha.request.json-template" -> """{"url": "${drsPath}"}"""
    ).asJava
  )

  private val marthaV3Config: Config = ConfigFactory.parseMap(
    Map(
      "martha.url" -> "https://martha-url/martha_v3",
      "martha.request.json-template" -> """{"url": "${drsPath}"}"""
    ).asJava
  )

  private lazy val fakeCredentials = NoCredentials.getInstance

  private lazy val httpClientBuilder = HttpClientBuilder.create()

  private def drsReadInterpreter(marthaResponse: MarthaResponse): IO[ReadableByteChannel] =
    throw new UnsupportedOperationException("Currently DrsResolverSpec doesn't need to use drs read interpreter.")

  private val mockFileSystemProviderForMarthaV2 = new MockDrsCloudNioFileSystemProvider(marthaV2Config, fakeCredentials, httpClientBuilder, drsReadInterpreter)
  private val drsPathBuilderForMarthaV2 = DrsPathBuilder(mockFileSystemProviderForMarthaV2, None)

  private val mockFileSystemProviderForMarthaV3 = new MockDrsCloudNioFileSystemProvider(marthaV3Config, fakeCredentials, httpClientBuilder, drsReadInterpreter)
  private val drsPathBuilderForMarthaV3 = DrsPathBuilder(mockFileSystemProviderForMarthaV3, None)


  behavior of "DrsResolver"

  it should "find GCS path" in {
    val drsPath = drsPathBuilderForMarthaV2.build(MockDrsPaths.drsPathResolvingGcsPath).get.asInstanceOf[DrsPath]

    DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync() should be (MockDrsPaths.gcsRelativePath)
  }

  it should "throw GcsUrlNotFoundException when DRS path doesn't resolve to at least one GCS url" in {
    val drsPath = drsPathBuilderForMarthaV2.build(MockDrsPaths.drsPathResolvingToNoGcsPath).get.asInstanceOf[DrsPath]

    the[RuntimeException] thrownBy {
      DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
    } should have message s"Error while resolving DRS path: ${drsPath.pathAsString}. Error: UrlNotFoundException: No gs url associated with given DRS path."
  }

  it should "throw Runtime Exception when Martha V2 can't find the given DRS path " in {
    val drsPath = drsPathBuilderForMarthaV2.build(MockDrsPaths.drsPathNotExistingInMartha).get.asInstanceOf[DrsPath]

    the[RuntimeException] thrownBy {
      DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
    } should have message s"Unexpected response resolving ${drsPath.pathAsString} through Martha url https://martha-url/martha_v2. Error: 404 Not Found."
  }

  it should "throw Runtime Exception when Martha V3 can't find the given DRS path " in {
    val drsPath = drsPathBuilderForMarthaV3.build(MockDrsPaths.drsPathNotExistingInMartha).get.asInstanceOf[DrsPath]

    the[RuntimeException] thrownBy {
      DrsResolver.getContainerRelativePath(drsPath).unsafeRunSync()
    } should have message s"Unexpected response resolving ${drsPath.pathAsString} through Martha url https://martha-url/martha_v3. Error: 404 Not Found."
  }
}
