package cromwell.engine.io.gcs

import com.google.cloud.storage.StorageException
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers


class GcsBatchFlowSpec extends AnyFlatSpecLike with Matchers {

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
    storageExceptions(objectReadOperationNames) map GcsBatchFlow.getReadForbiddenBucket map { _.get } shouldBe Set(UnreadableBucketName)

    // A sampling of write operations, not an exhaustive list.
    val objectWriteOperationNames = Set(
      "compose",
      "delete",
      "insert"
    )
    storageExceptions(objectWriteOperationNames) flatMap GcsBatchFlow.getReadForbiddenBucket shouldBe Set.empty

    Set(new RuntimeException("random exception")) flatMap GcsBatchFlow.getReadForbiddenBucket shouldBe Set.empty
  }
}
