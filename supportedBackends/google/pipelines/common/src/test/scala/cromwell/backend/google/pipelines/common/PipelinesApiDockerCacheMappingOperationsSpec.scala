package cromwell.backend.google.pipelines.common

import cats.effect.IO
import com.google.cloud.storage.{Blob, BlobId, Storage}
import common.assertion.CromwellTimeoutSpec
import cromwell.core.logging.JobLogger
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

  private val pipelinesApiDockerCacheMappingOperationsMock = new PipelinesApiDockerCacheMappingOperations {}

  it should "successfully parse docker image cache manifest JSON file as instance of Map[String, String]" in {
    val expectedManifest = DockerImageCacheManifest(
      manifestFormatVersion = 2,
      dockerImageCacheMap =
        Map(
          "project1/dockerImage1" ->
            DockerImageCacheEntry(
              "sha256:c73e11b1a7854f31ff12c607738bef7b0560a880dc9d0445dde084acc0e9da09",
              "projects/some-google-project/global/images/dockerCacheDiskForDockerImage1"
            ),
          "project2/dockerImage2" ->
            DockerImageCacheEntry(
              "sha256:7eb386481d87e41ebddceb948f25379bb339784df5247a3dfbdea2ac101b10c0",
              "projects/another-google-project/global/images/dockerCacheDiskForDockerImage2"
            )
        )
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

    val readFileFromGcsPrivateMethod = PrivateMethod[IO[DockerImageCacheManifest]]('readDockerImageCacheManifestFileFromGCS)
    val parsedJsonAsManifestIO = pipelinesApiDockerCacheMappingOperationsMock invokePrivate readFileFromGcsPrivateMethod(mockGcsClient, testJsonGcsPath)
    val parsedJsonAsManifest = parsedJsonAsManifestIO.unsafeRunSync()

    parsedJsonAsManifest.equals(expectedManifest) shouldBe true
  }

  it should "not use docker image cache if digest of requested image and cached images don't match" in {
    val mockJobLogger = mock[JobLogger]
    val testDockerImageName = "test_madeup_docker_image"
    val testDockerImageDigest = "fake_docker_image_digest"
    val testDiskImageName = "fake_disk_image_name"
    val dockerImageToCacheDiskImageMapping = Map(testDockerImageName -> DockerImageCacheEntry(testDockerImageDigest, testDiskImageName))
    val dockerImageCacheDiskOpt =
      pipelinesApiDockerCacheMappingOperationsMock.getDockerCacheDiskImageForAJob(
        dockerImageToCacheDiskImageMappingOpt = Option(dockerImageToCacheDiskImageMapping),
        dockerImageAsSpecifiedByUser = testDockerImageName,
        dockerImageWithDigest = s"$testDockerImageName@not_matching_digest",
        jobLogger = mockJobLogger
      )

    dockerImageCacheDiskOpt shouldBe None
  }

  it should "use docker image cache if digest of requested image and cached image match" in {
    val mockJobLogger = mock[JobLogger]
    val testDockerImageName = "test_madeup_docker_image"
    val testDockerImageDigest = "fake_docker_image_digest"
    val expectedDiskImageName = "fake_disk_image_name"
    val dockerImageToCacheDiskImageMapping = Map(testDockerImageName -> DockerImageCacheEntry(testDockerImageDigest, expectedDiskImageName))
    val dockerImageCacheDiskOpt =
      pipelinesApiDockerCacheMappingOperationsMock.getDockerCacheDiskImageForAJob(
        dockerImageToCacheDiskImageMappingOpt = Option(dockerImageToCacheDiskImageMapping),
        dockerImageAsSpecifiedByUser = testDockerImageName,
        dockerImageWithDigest = s"$testDockerImageName@$testDockerImageDigest",
        jobLogger = mockJobLogger
      )

    dockerImageCacheDiskOpt shouldBe Some(expectedDiskImageName)
  }

  it should "not use docker image cache if requested docker image is not in cache" in {
    val mockJobLogger = mock[JobLogger]
    val testCachedDockerImageName = "test_madeup_docker_image"
    val testDockerImageNameSpecifiedByUser = "test_non_cached_docker_image"
    val testDockerImageDigest = "fake_docker_image_digest"
    val testDiskImageName = "fake_disk_image_name"
    val dockerImageToCacheDiskImageMapping = Map(testCachedDockerImageName -> DockerImageCacheEntry(testDockerImageDigest, testDiskImageName))
    val dockerImageCacheDiskOpt =
      pipelinesApiDockerCacheMappingOperationsMock.getDockerCacheDiskImageForAJob(
        dockerImageToCacheDiskImageMappingOpt = Option(dockerImageToCacheDiskImageMapping),
        dockerImageAsSpecifiedByUser = testDockerImageNameSpecifiedByUser,
        dockerImageWithDigest = s"$testDockerImageNameSpecifiedByUser@$testDockerImageDigest",
        jobLogger = mockJobLogger
      )

    dockerImageCacheDiskOpt shouldBe None
  }

}
