package cromwell.engine.io

import java.util.UUID

import akka.stream.ActorMaterializer
import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import com.typesafe.config.ConfigFactory
import cromwell.core.Tags.IntegrationTest
import cromwell.core.io._
import cromwell.core.{TestKitSuite, WorkflowOptions}
import cromwell.filesystems.gcs.batch._
import cromwell.filesystems.gcs.{GcsPathBuilder, GcsPathBuilderFactory}
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.concurrent.{Await, ExecutionContext}
import scala.language.postfixOps

class IoActorProxyGcsBatchSpec extends TestKitSuite with FlatSpecLike with Matchers with ImplicitSender with Eventually {
  behavior of "IoActor [GCS Batch]"

  implicit val actorSystem = system
  implicit val ec: ExecutionContext = system.dispatcher
  implicit val materializer = ActorMaterializer()
  
  val instanceConfig = ConfigFactory.parseString("auth = \"application-default\"")

  override def afterAll() = {
    materializer.shutdown()
    src.delete(swallowIOExceptions = true)
    dst.delete(swallowIOExceptions = true)
    srcRegional.delete(swallowIOExceptions = true)
    dstMultiRegional.delete(swallowIOExceptions = true)
    super.afterAll()
  }
  
  lazy val gcsPathBuilder = GcsPathBuilderFactory(ConfigFactory.load(), instanceConfig)
  lazy val pathBuilder: GcsPathBuilder = Await.result(gcsPathBuilder.withOptions(WorkflowOptions.empty), 1 second)

  lazy val randomUUID = UUID.randomUUID().toString

  lazy val directory = pathBuilder.build(s"gs://cloud-cromwell-dev/unit-test").get
  lazy val src = pathBuilder.build(s"gs://cloud-cromwell-dev/unit-test/$randomUUID/testFile.txt").get
  lazy val dst = pathBuilder.build(s"gs://cloud-cromwell-dev/unit-test/$randomUUID/testFile-copy.txt").get
  lazy val srcRegional = pathBuilder.build(s"gs://cloud-cromwell-dev-regional/unit-test/$randomUUID/testRegional.txt").get
  lazy val dstMultiRegional = pathBuilder.build(s"gs://cloud-cromwell-dev/unit-test/$randomUUID/testFileRegional-copy.txt").get
  
  override def beforeAll() = {
    // Write commands can't be batched, so for the sake of this test, just create a file in GCS synchronously here
    src.write("hello")
    srcRegional.write("hello")
    super.beforeAll()
  }
  
  it should "batch queries" taggedAs IntegrationTest in {
    val testActor = TestActorRef(new IoActor(10, 10, 10, None, TestProbe().ref))

    val copyCommand = GcsBatchCopyCommand(src, dst, overwrite = false)
    val sizeCommand = GcsBatchSizeCommand(src)
    val hashCommand = GcsBatchCrc32Command(src)
    // Should return true
    val isDirectoryCommand = GcsBatchIsDirectoryCommand(directory)
    // Should return false
    val isDirectoryCommand2 = GcsBatchIsDirectoryCommand(src)
    
    val deleteSrcCommand = GcsBatchDeleteCommand(src, swallowIOExceptions = false)
    val deleteDstCommand = GcsBatchDeleteCommand(dst, swallowIOExceptions = false)

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

  it should "copy files across GCS storage classes" taggedAs IntegrationTest in {
    val testActor = TestActorRef(new IoActor(10, 10, 10, None, TestProbe().ref))

    val copyCommand = GcsBatchCopyCommand(srcRegional, dstMultiRegional, overwrite = false)
    
    testActor ! copyCommand

    expectMsgClass(30 seconds, classOf[IoSuccess[_]])

    dstMultiRegional.exists shouldBe true
  }
}
