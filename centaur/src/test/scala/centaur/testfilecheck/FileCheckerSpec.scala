package centaur.testfilecheck

import java.util

import com.google.api.gax.paging.Page
import com.google.cloud.storage.Storage.BlobListOption
import com.google.cloud.storage.{Blob, Storage}
import org.mockito.Mockito._
import org.scalatest.FlatSpec
import org.specs2.matcher.ShouldMatchers
import org.specs2.mock.Mockito
import software.amazon.awssdk.services.s3.S3Client
import software.amazon.awssdk.services.s3.model.{ListObjectsRequest, ListObjectsResponse, S3Object}

class FileCheckerSpec extends FlatSpec with ShouldMatchers with Mockito {

  import centaur.test.ObjectCounterInstances._

  val s3PrefixRegex = "^s3:\\/\\/.*"
  val gsPrefixRegex = "^gs:\\/\\/.*"

  val amazonS3mock = mock[S3Client]
  val testPath = "s3://my-cool-bucket/path/to/file"
  val bucketName = "my-cool-bucket"
  val dirName = "path/to/file"
  val wrongBucketPrefix = "s3Bucket://my-not-so-cool-bucket/somelogs/empty"
  val EmptyTestPath = ""
  val gsPathType = "gs://"
  val testGsPath = "gs://my-cool-bucket/path/to/file"
  val objResponse = ListObjectsResponse.builder()
    .contents(util.Arrays.asList(S3Object.builder()
      .build()))
    .build()
  val objRequest = ListObjectsRequest.builder().bucket(bucketName).prefix(dirName).build()
  val awsS3Path = awsS3ObjectCounter.parsePath(s3PrefixRegex)(testPath)
  val gsPath = gcsObjectCounter.parsePath(gsPrefixRegex)(testGsPath)


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
    when(pageMock.iterateAll).thenReturns(blobbies)
    when(gcsMock.list(bucketName, BlobListOption.prefix(dirName))).thenReturn(pageMock)
    val actualObjectCounts = gcsObjectCounter.countObjectsAtPath(gcsMock)(gsPath)
    val expectedObjectCounts = 1
    assert(expectedObjectCounts == actualObjectCounts)
  }

  "countObjectsAtPath GS" should "return 0 for empty list" in {

    val gcsMock = mock[Storage]
    val pageMock = mock[Page[Blob]]
    val blobbies = new util.ArrayList[Blob]()
    when(pageMock.iterateAll).thenReturns(blobbies)
    when(gcsMock.list(bucketName, BlobListOption.prefix(dirName))).thenReturn(pageMock)
    val actualObjectCounts = gcsObjectCounter.countObjectsAtPath(gcsMock)(gsPath)
    val expectedObjectCounts = 0
    assert(expectedObjectCounts == actualObjectCounts)
  }
}