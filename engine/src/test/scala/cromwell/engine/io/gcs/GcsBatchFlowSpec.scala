package cromwell.engine.io.gcs

import akka.testkit.TestProbe
import com.google.cloud.storage.StorageException
import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem
import common.assertion.CromwellTimeoutSpec
import cromwell.core.TestKitSuite
import cromwell.filesystems.gcs.GcsPath
import cromwell.filesystems.gcs.batch.GcsBatchCrc32Command
import org.specs2.mock.Mockito
import org.scalatest.PrivateMethodTester
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.{ExecutionContextExecutor, Future}


class GcsBatchFlowSpec extends TestKitSuite with AnyFlatSpecLike with CromwellTimeoutSpec with Matchers with PrivateMethodTester with Mockito {

  "GcsBatchFlow" should "know what read forbidden bucket failures look like" in {
    val ErrorTemplate = "foo@bar.iam.gserviceaccount.com does not have storage.objects.%s access to %s/three_step/f0000000-baaa-f000-baaa-f00000000000/call-foo/foo.log"
    val UnreadableBucketName = "unreadable-bucket"

    val objectReadOperationNames = Set(
      "list",
      "get",
      "copy"
    )

    def storageExceptions(names: Set[String]): Set[StorageException] = {
      val failureCodes = List(403, 404)
      for {
        op <- names
        code <- failureCodes
      } yield new StorageException(code, String.format(ErrorTemplate, op, UnreadableBucketName))
    }
    // Can't flatMap this since any Nones would just be squashed out.
    storageExceptions(objectReadOperationNames) map (_.getMessage) map GcsBatchFlow.getReadForbiddenBucket map { _.get } shouldBe Set(UnreadableBucketName)

    // A sampling of write operations, not an exhaustive list.
    val objectWriteOperationNames = Set(
      "compose",
      "delete",
      "insert"
    )
    storageExceptions(objectWriteOperationNames) map (_.getMessage) flatMap GcsBatchFlow.getReadForbiddenBucket shouldBe Set.empty

    Set(new RuntimeException("random exception")) map (_.getMessage) flatMap GcsBatchFlow.getReadForbiddenBucket shouldBe Set.empty
  }

  "GcsBatchFlow" should "not throw unhandled exception and kill the thread when trying to recover from unretryable exception with null error message" in {
    implicit val ec: ExecutionContextExecutor = system.dispatcher
    val gcsBatchFlow = new GcsBatchFlow(1, system.scheduler, null, "testAppName")

    val mockGcsPath = GcsPath(CloudStorageFileSystem.forBucket("bucket").getPath("test"), any[com.google.api.services.storage.Storage], any[com.google.cloud.storage.Storage], anyString)
    val gcsBatchCommandContext = GcsBatchCommandContext(GcsBatchCrc32Command(mockGcsPath), TestProbe().ref)
    val recoverCommandPrivateMethod = PrivateMethod[PartialFunction[Throwable, Future[GcsBatchResponse[_]]]]('recoverCommand)
    val partialFuncAcceptingThrowable = gcsBatchFlow invokePrivate recoverCommandPrivateMethod(gcsBatchCommandContext)

    val futureRes = partialFuncAcceptingThrowable(new NullPointerException(null)) // no unhandled exceptions should be thrown here
    futureRes.isCompleted shouldBe true
  }
}
