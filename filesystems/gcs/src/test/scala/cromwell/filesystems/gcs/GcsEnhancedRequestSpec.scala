package cromwell.filesystems.gcs

import java.io.FileNotFoundException

import com.google.api.client.googleapis.json.{GoogleJsonError, GoogleJsonResponseException}
import com.google.api.client.http.{HttpHeaders, HttpResponseException}
import com.google.cloud.storage.StorageException
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import cromwell.filesystems.gcs.RequesterPaysErrors._
import org.scalamock.scalatest.MockFactory
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.specs2.mock.Mockito

class GcsEnhancedRequestSpec extends AnyFlatSpec with Matchers with Mockito with MockFactory {
  behavior of "GcsEnhancedRequest"

  val path = GcsPath(CloudStorageFileSystem.forBucket("bucket").getPath("test"), any[com.google.api.services.storage.Storage], any[com.google.cloud.storage.Storage], anyString)
  val requesterPaysException = new StorageException(BucketIsRequesterPaysErrorCode, BucketIsRequesterPaysErrorMessage)
  
  it should "attempt first without project, and not retry if the requests succeeds" in {
    val testFunction = mockFunction[Boolean, String]
    testFunction.expects(false).returns("hello").once()
    GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync() shouldBe "hello"
  }

  it should "retry requests with project if the error matches and succeed" in {
    val testFunction = mockFunction[Boolean, String]

    // First time, throw a requester pays exception
    testFunction.expects(false).throws(requesterPaysException)
    // We expect it to be called a second time with withProject = true this time
    testFunction.expects(true).returns("hello")
    GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync() shouldBe "hello"
  }

  it should "retry requests with project if the error matches and fail" in {
    val testFunction = mockFunction[Boolean, String]

    // First time, throw a requester pays exception
    testFunction.expects(false).throws(requesterPaysException)
    // We expect it to be called a second time with withProject = true this time, and fail for another reason
    testFunction.expects(true).throws(new RuntimeException("it really doesn't work"))
    a[RuntimeException] should be thrownBy GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync()
  }

  it should "not retry requests if the error does not match" in {
    val testFunction = mockFunction[Boolean, String]

    // Throw an unrelated exception, should only be called once
    testFunction.expects(false).throws(new RuntimeException("it really doesn't work")).once()
    a[RuntimeException] should be thrownBy GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync()
  }

  it should "re throw StorageException 404 to NoSuchFileException" in {
    val testFunction = mockFunction[Boolean, String]

    // Throw an unrelated exception, should only be called once
    testFunction.expects(false).throws(new StorageException(404, "gs://does/not/exist")).once()
    a[FileNotFoundException] should be thrownBy GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync()
  }

  it should "re throw GoogleJsonResponseException 404 to NoSuchFileException" in {
    val testFunction = mockFunction[Boolean, String]

    val builder = new HttpResponseException.Builder(404, "NotFound", new HttpHeaders)
    val error = new GoogleJsonError()
    error.setCode(404)
    
    // Throw an unrelated exception, should only be called once
    testFunction.expects(false).throws(new GoogleJsonResponseException(builder, error)).once()
    a[FileNotFoundException] should be thrownBy GcsEnhancedRequest.recoverFromProjectNotProvided(path, testFunction).unsafeRunSync()
  }
}
