package cromwell.filesystems.gcs

import java.io.FileNotFoundException

import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonResponseException}
import com.google.api.client.http.{HttpHeaders, HttpResponseException}
import com.google.cloud.storage.StorageException
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import common.assertion.CromwellTimeoutSpec
import cromwell.filesystems.gcs.RequesterPaysErrors._
import org.mockito.ArgumentMatchers._
import org.mockito.Mockito._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import common.mock.MockSugar

class GcsEnhancedRequestSpec extends AnyFlatSpec with CromwellTimeoutSpec with Matchers with MockSugar {
  behavior of "GcsEnhancedRequest"

  private val path =
    GcsPath(
      CloudStorageFileSystem.forBucket("bucket").getPath("test"),
      mock[com.google.api.services.storage.Storage],
      mock[com.google.cloud.storage.Storage],
      "GcsEnhancedRequest-project",
    )
  val requesterPaysException = new StorageException(BucketIsRequesterPaysErrorCode, "Bucket is a requester pays bucket but no user project provided.")

  it should "attempt first without project, and not retry if the requests succeeds" in {
    val testFunction = mock[Boolean => String]
    when(testFunction.apply(false)).thenReturn("hello")
    GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync() shouldBe "hello"
    verify(testFunction).apply(false)
    verify(testFunction).apply(anyBoolean)
  }

  it should "retry requests with project if the error matches and succeed" in {
    val testFunction = mock[Boolean => String]

    // First time, throw a requester pays exception
    when(testFunction.apply(false)).thenThrow(requesterPaysException)
    // We expect it to be called a second time with withProject = true this time
    when(testFunction.apply(true)).thenReturn("hello")
    GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync() shouldBe "hello"
    verify(testFunction).apply(false)
    verify(testFunction).apply(true)
    verify(testFunction, times(2)).apply(anyBoolean)
  }

  it should "retry requests with project if the error matches and fail" in {
    val testFunction = mock[Boolean => String]

    // First time, throw a requester pays exception
    when(testFunction.apply(false)).thenThrow(requesterPaysException)
    // We expect it to be called a second time with withProject = true this time, and fail for another reason
    when(testFunction.apply(true)).thenThrow(new RuntimeException("it really doesn't work"))
    a[RuntimeException] should be thrownBy GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync()
  }

  it should "not retry requests if the error does not match" in {
    val testFunction = mock[Boolean => String]

    // Throw an unrelated exception, should only be called once
    when(testFunction.apply(false)).thenThrow(new RuntimeException("it really doesn't work"))
    a[RuntimeException] should be thrownBy GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync()
    verify(testFunction).apply(false)
    verify(testFunction).apply(anyBoolean)
  }

  it should "re throw StorageException 404 to NoSuchFileException" in {
    val testFunction = mock[Boolean => String]

    // Throw an unrelated exception, should only be called once
    when(testFunction.apply(false)).thenThrow(new StorageException(404, "gs://does/not/exist"))
    a[FileNotFoundException] should be thrownBy GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync()
    verify(testFunction).apply(false)
    verify(testFunction).apply(anyBoolean)
  }

  it should "re throw GoogleJsonResponseException 404 to NoSuchFileException" in {
    val testFunction = mock[Boolean => String]

    val builder = new HttpResponseException.Builder(404, "NotFound", new HttpHeaders)
    val error = new GoogleJsonError()
    error.setCode(404)

    // Throw an unrelated exception, should only be called once
    when(testFunction.apply(false)).thenAnswer(_ => throw new GoogleJsonResponseException(builder, error))
    a[FileNotFoundException] should be thrownBy GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync()
    verify(testFunction).apply(false)
    verify(testFunction).apply(anyBoolean)
  }
}
