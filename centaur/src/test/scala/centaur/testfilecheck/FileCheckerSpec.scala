package centaur.testfilecheck

import java.util
import com.google.api.gax.paging.Page
import com.google.cloud.storage.Storage.BlobListOption
import com.google.cloud.storage.{Blob, Storage}
import common.assertion.CromwellTimeoutSpec
import org.mockito.Mockito._
import common.mock.MockSugar
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.model.{ListObjectsRequest, ListObjectsResponse, S3Object}
import org.scalatest.flatspec.AnyFlatSpec


class FileCheckerSpec extends AnyFlatSpec with CromwellTimeoutSpec with MockSugar {

  import centaur.test.ObjectCounterInstances._

  private val s3PrefixRegex = "^s3:\\/\\/.*"
  private val gsPrefixRegex = "^gs:\\/\\/.*"

  private val amazonS3mock = mock[S3Client]
  private val testPath = "s3://my-cool-bucket/path/to/file"
  private val bucketName = "my-cool-bucket"
  private val dirName = "path/to/file"
  private val wrongBucketPrefix = "s3Bucket://my-not-so-cool-bucket/somelogs/empty"
  private val EmptyTestPath = ""
  private val testGsPath = "gs://my-cool-bucket/path/to/file"
  private val objResponse = ListObjectsResponse.builder()
    .contents(util.Arrays.asList(S3Object.builder()
      .build()))
    .build()
  private val objRequest = ListObjectsRequest.builder().bucket(bucketName).prefix(dirName).build()
  private val awsS3Path = awsS3ObjectCounter.parsePath(s3PrefixRegex)(testPath)
  private val gsPath = gcsObjectCounter.parsePath(gsPrefixRegex)(testGsPath)


  "parsePath" should "return a bucket and directories" in {
    assert(awsS3Path.bucket == bucketName)
    assert(awsS3Path.directory == dirName)
  }

  "parsePath" should "throw Exception for wrong path" in {
    assertThrows[centaur.test.IllegalPathException] {awsS3ObjectCounter.parsePath(s3PrefixRegex)(wrongBucketPrefix)}
    assertThrows[centaur.test.IllegalPathException] {awsS3ObjectCounter.parsePath(s3PrefixRegex)(testGsPath)}
    assertThrows[centaur.test.IllegalPathException] {awsS3ObjectCounter.parsePath(s3PrefixRegex)(EmptyTestPath)}
  }

  "countObjectAtPath" should "should return 1 if the file exist" in {
    when(amazonS3mock.listObjects(objRequest)).thenReturn(objResponse)
    val actualObjectCounts = awsS3ObjectCounter.countObjectsAtPath(amazonS3mock)(awsS3Path)
    val expectedObjectCounts = 1
    assert(expectedObjectCounts == actualObjectCounts)
  }

  "countObjectsAtPath" should "return 0 for empty list" in {

    val objResponse = ListObjectsResponse.builder().build()
    when(amazonS3mock.listObjects(objRequest)).thenReturn(objResponse)
    val actualObjectCounts = awsS3ObjectCounter.countObjectsAtPath(amazonS3mock)(awsS3Path)
    val expectedObjectCounts = 0
    assert(expectedObjectCounts == actualObjectCounts)
  }

  "gcsObjectCounter.parsePath" should "return a bucket and directories" in {
    assert(gsPath.bucket == bucketName)
    assert(gsPath.directory == dirName)
  }

  "gcsObjectCounter.parsePath" should "return 1 if the file exist" in {

    val gcsMock = mock[Storage]
    val pageMock = mock[Page[Blob]]
    val blob = null
    val blobbies = new util.ArrayList[Blob]()
    blobbies.add(blob)
    when(pageMock.iterateAll).thenReturn(blobbies)
    when(gcsMock.list(bucketName, BlobListOption.prefix(dirName))).thenReturn(pageMock)
    val actualObjectCounts = gcsObjectCounter.countObjectsAtPath(gcsMock)(gsPath)
    val expectedObjectCounts = 1
    assert(expectedObjectCounts == actualObjectCounts)
  }

  "countObjectsAtPath GS" should "return 0 for empty list" in {

    val gcsMock = mock[Storage]
    val pageMock = mock[Page[Blob]]
    val blobbies = new util.ArrayList[Blob]()
    when(pageMock.iterateAll).thenReturn(blobbies)
    when(gcsMock.list(bucketName, BlobListOption.prefix(dirName))).thenReturn(pageMock)
    val actualObjectCounts = gcsObjectCounter.countObjectsAtPath(gcsMock)(gsPath)
    val expectedObjectCounts = 0
    assert(expectedObjectCounts == actualObjectCounts)
  }
}
