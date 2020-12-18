package cromwell.backend.google.pipelines.common

import cats.effect.IO
import com.google.cloud.storage.{Blob, BlobId, Storage}
import common.assertion.CromwellTimeoutSpec
import cromwell.filesystems.gcs.GcsPathBuilder.ValidFullGcsPath
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.mockito.Mockito._
import org.scalatest.PrivateMethodTester
import org.specs2.mock.Mockito

import scala.io.Source

class PipelinesApiDockerCacheMappingOperationsSpec
  extends AnyFlatSpecLike
    with CromwellTimeoutSpec
    with Matchers
    with Mockito
    with PrivateMethodTester {

  it should "successfully parse docker image cache manifest JSON file as instance of Map[String, String]" in {
    val pipelinesApiDockerCacheMappingOperationsMock = new PipelinesApiDockerCacheMappingOperations {}
    val expectedMap = Map(
      "project1/dockerImage1" -> "projects/some-google-project/global/images/dockerCacheDiskForDockerImage1",
      "project2/dockerImage2" -> "projects/another-google-project/global/images/dockerCacheDiskForDockerImage2"
    )

    val testJsonFileName = "docker-image-cache-manifest.json"
    val testJsonGcsPath = ValidFullGcsPath("test-bucket", s"/$testJsonFileName")
    val mockBlobId = BlobId.of(testJsonGcsPath.bucket, testJsonGcsPath.path.substring(1))

    val mockJsonBlob = {
      val mockBlob = mock[Blob]
      val testJsonAsByteArray = Source.fromInputStream(
        Thread.currentThread.getContextClassLoader.getResourceAsStream(testJsonFileName)
      ).mkString.getBytes

      when(mockBlob.getContent()).thenReturn(testJsonAsByteArray)
      mockBlob
    }

    val mockGcsClient = {
      val mockClient = mock[Storage]
      when(mockClient.get(mockBlobId)).thenReturn(mockJsonBlob)
      mockClient
    }

    val readFileFromGcsPrivateMethod = PrivateMethod[IO[Map[String, String]]]('readDockerImageCacheManifestFileFromGCS)
    val parsedJsonAsMapIO = pipelinesApiDockerCacheMappingOperationsMock invokePrivate readFileFromGcsPrivateMethod(mockGcsClient, testJsonGcsPath)
    val parsedJsonAsMap = parsedJsonAsMapIO.unsafeRunSync()

    parsedJsonAsMap.equals(expectedMap) shouldBe true
  }

}
