package cromwell.engine.io

import java.util.UUID

import akka.stream.ActorMaterializer
import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import com.typesafe.config.ConfigFactory
import cromwell.core.Tags.IntegrationTest
import cromwell.core.io._
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.filesystems.gcs.batch._
import cromwell.filesystems.gcs.{GcsPath, GcsPathBuilder, GcsPathBuilderFactory}
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext}
import scala.language.postfixOps

class IoActorProxyGcsBatchSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with ImplicitSender with Eventually {
  behavior of "IoActor [GCS Batch]"

  implicit val ec: ExecutionContext = system.dispatcher
  implicit val materializer: ActorMaterializer = ActorMaterializer()

  private val instanceConfig = ConfigFactory.parseString("auth = \"integration-test\"")

  override def afterAll(): Unit = {
    materializer.shutdown()
    src.delete(swallowIOExceptions = true)
    srcRequesterPays.delete(swallowIOExceptions = true)
    dst.delete(swallowIOExceptions = true)
    dstRequesterPays.delete(swallowIOExceptions = true)
    srcRegional.delete(swallowIOExceptions = true)
    dstMultiRegional.delete(swallowIOExceptions = true)
    super.afterAll()
  }

  private lazy val gcsPathBuilder = GcsPathBuilderFactory(ConfigFactory.load(), instanceConfig)
  private lazy val pathBuilder: GcsPathBuilder =
    Await.result(gcsPathBuilder.withOptions(WorkflowOptions.empty), 30.seconds)

  private lazy val randomUUID = UUID.randomUUID().toString

  private lazy val directory = pathBuilder.build(s"gs://cloud-cromwell-dev/unit-test").get
  private lazy val src = pathBuilder.build(s"gs://cloud-cromwell-dev/unit-test/$randomUUID/testFile.txt").get
  private lazy val dst = pathBuilder.build(s"gs://cloud-cromwell-dev/unit-test/$randomUUID/testFile-copy.txt").get

  private lazy val directoryRequesterPays =
    pathBuilder.build(s"gs://cromwell_bucket_with_requester_pays/unit-test").get
  private lazy val srcRequesterPays =
    pathBuilder.build(s"gs://cromwell_bucket_with_requester_pays/unit-test/$randomUUID/testFile.txt").get
  private lazy val dstRequesterPays =
    pathBuilder.build(s"gs://cromwell_bucket_with_requester_pays/unit-test/$randomUUID/testFile-copy.txt").get

  private lazy val srcRegional =
    pathBuilder.build(s"gs://cloud-cromwell-dev-regional/unit-test/$randomUUID/testRegional.txt").get
  private lazy val dstMultiRegional =
    pathBuilder.build(s"gs://cloud-cromwell-dev/unit-test/$randomUUID/testFileRegional-copy.txt").get

  override def beforeAll(): Unit = {
    // Write commands can't be batched, so for the sake of this test, just create a file in GCS synchronously here
    src.write("hello")
    srcRequesterPays.write("hello")
    srcRegional.write("hello")
    super.beforeAll()
  }

  private def testWith(src: GcsPath,
                       dst: GcsPath,
                       directory: GcsPath,
                       testActorName: String,
                       serviceRegistryActorName: String) = {
    val testActor = TestActorRef(
      factory = new IoActor(10, 10, 10, None, TestProbe(serviceRegistryActorName).ref, "cromwell test"),
      name = testActorName,
    )

    val copyCommand = GcsBatchCopyCommand.forPaths(src, dst).get
    val sizeCommand = GcsBatchSizeCommand.forPath(src).get
    val hashCommand = GcsBatchCrc32Command.forPath(src).get
    // Should return true
    val isDirectoryCommand = GcsBatchIsDirectoryCommand.forPath(directory).get
    // Should return false
    val isDirectoryCommand2 = GcsBatchIsDirectoryCommand.forPath(src).get

    val deleteSrcCommand = GcsBatchDeleteCommand.forPath(src, swallowIOExceptions = false).get
    val deleteDstCommand = GcsBatchDeleteCommand.forPath(dst, swallowIOExceptions = false).get

    testActor ! copyCommand
    testActor ! sizeCommand
    testActor ! hashCommand
    testActor ! isDirectoryCommand
    testActor ! isDirectoryCommand2

    val received1 = receiveN(5, 10 seconds)

    received1.size shouldBe 5
    received1 forall { _.isInstanceOf[IoSuccess[_]] } shouldBe true

    received1 collect {
      case IoSuccess(_: GcsBatchSizeCommand, fileSize: Long) => fileSize shouldBe 5
    }

    received1 collect {
      case IoSuccess(_: GcsBatchCrc32Command, hash: String) => hash shouldBe "mnG7TA=="
    }

    received1 collect {
      case IoSuccess(command: GcsBatchIsDirectoryCommand, isDirectory: Boolean) if command.file.pathAsString == directory.pathAsString => isDirectory shouldBe true
    }

    received1 collect {
      case IoSuccess(command: GcsBatchIsDirectoryCommand, isDirectory: Boolean) if command.file.pathAsString == src.pathAsString => isDirectory shouldBe false
    }

    testActor ! deleteSrcCommand
    testActor ! deleteDstCommand

    val received2 = receiveN(2, 10 seconds)

    received2.size shouldBe 2
    received2 forall { _.isInstanceOf[IoSuccess[_]] } shouldBe true

    src.exists shouldBe false
    dst.exists shouldBe false
  }

  it should "batch queries" taggedAs IntegrationTest in {
    testWith(
      src = src,
      dst = dst,
      directory = directory,
      testActorName = "testActor-batch",
      serviceRegistryActorName = "serviceRegistryActor-batch",
    )
  }

  it should "batch queries on requester pays buckets" taggedAs IntegrationTest in {
    testWith(
      src = srcRequesterPays,
      dst = dstRequesterPays,
      directory = directoryRequesterPays,
      testActorName = "testActor-batch-rp",
      serviceRegistryActorName = "serviceRegistryActor-batch-rp",
    )
  }

  it should "copy files across GCS storage classes" taggedAs IntegrationTest in {
    val testActor = TestActorRef(
      factory = new IoActor(10, 10, 10, None, TestProbe("serviceRegistryActor").ref, "cromwell test"),
      name = "testActor",
    )

    val copyCommand = GcsBatchCopyCommand.forPaths(srcRegional, dstMultiRegional).get

    testActor ! copyCommand

    expectMsgClass(30 seconds, classOf[IoSuccess[_]])

    dstMultiRegional.exists shouldBe true
  }
}
