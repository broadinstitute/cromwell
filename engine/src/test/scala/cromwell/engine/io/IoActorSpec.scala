package cromwell.engine.io

import java.io.IOException
import java.net.{SocketException, SocketTimeoutException}

import akka.stream.ActorMaterializer
import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import better.files.File.OpenOptions
import com.google.cloud.storage.StorageException
import cromwell.core.TestKitSuite
import cromwell.core.io.DefaultIoCommand._
import cromwell.core.io.IoContentAsStringCommand.IoReadOptions
import cromwell.core.io._
import cromwell.core.path.{DefaultPathBuilder, Path}
import cromwell.engine.io.gcs.GcsBatchFlow.BatchFailedException
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps

class IoActorSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with ImplicitSender {
  behavior of "IoActor"
  
  implicit val ec: ExecutionContext = system.dispatcher
  implicit val materializer: ActorMaterializer = ActorMaterializer()
  
  override def afterAll(): Unit = {
    materializer.shutdown()
    super.afterAll()
  }
  
  it should "copy a file" in {
    val testActor = TestActorRef(
      factory = new IoActor(1, 10, 10, None, TestProbe("serviceRegistryActorCopy").ref, "cromwell test"),
      name = "testActorCopy",
    )
    
    val src = DefaultPathBuilder.createTempFile()
    val dst: Path = src.parent.resolve(src.name + "-dst")
    
    val copyCommand = DefaultIoCopyCommand(src, dst)
    
    testActor ! copyCommand
    expectMsgPF(5 seconds) {
      case response: IoSuccess[_] => response.command.isInstanceOf[IoCopyCommand] shouldBe true
      case response: IoFailure[_] => fail("Expected an IoSuccess", response.failure)
    }
    
    dst.toFile should exist
    src.delete()
    dst.delete()
  }

  it should "write to a file" in {
    val testActor = TestActorRef(
      factory = new IoActor(1, 10, 10, None, TestProbe("serviceRegistryActorWrite").ref, "cromwell test"),
      name = "testActorWrite",
    )

    val src = DefaultPathBuilder.createTempFile()

    val writeCommand = DefaultIoWriteCommand(src, "hello", OpenOptions.default, compressPayload = false)

    testActor ! writeCommand
    expectMsgPF(5 seconds) {
      case response: IoSuccess[_] => response.command.isInstanceOf[IoWriteCommand] shouldBe true
      case response: IoFailure[_] => fail("Expected an IoSuccess", response.failure)
    }

    src.contentAsString shouldBe "hello"
    src.delete()
  }

  it should "delete a file" in {
    val testActor = TestActorRef(
      factory = new IoActor(1, 10, 10, None, TestProbe("serviceRegistryActorDelete").ref, "cromwell test"),
      name = "testActorDelete",
    )

    val src = DefaultPathBuilder.createTempFile()

    val deleteCommand = DefaultIoDeleteCommand(src, swallowIOExceptions = false)

    testActor ! deleteCommand
    expectMsgPF(5 seconds) {
      case response: IoSuccess[_] => response.command.isInstanceOf[IoDeleteCommand] shouldBe true
      case response: IoFailure[_] => fail("Expected an IoSuccess", response.failure)
    }

    src.toFile shouldNot exist
  }

  it should "read a file" in {
    val testActor = TestActorRef(
      factory = new IoActor(1, 10, 10, None, TestProbe("serviceRegistryActorRead").ref, "cromwell test"),
      name = "testActorRead",
    )

    val src = DefaultPathBuilder.createTempFile()
    src.write("hello")

    val readCommand = DefaultIoContentAsStringCommand(src, IoReadOptions(None, failOnOverflow = false))

    testActor ! readCommand
    expectMsgPF(5 seconds) {
      case response: IoSuccess[_] => 
        response.command.isInstanceOf[IoContentAsStringCommand] shouldBe true
        response.result.asInstanceOf[String] shouldBe "hello"
      case response: IoFailure[_] => fail("Expected an IoSuccess", response.failure)
    }

    src.delete()
  }

  it should "read only the first bytes of file" in {
    val testActor = TestActorRef(
      factory = new IoActor(1, 10, 10, None, TestProbe("serviceRegistryActorFirstBytes").ref, "cromwell test"),
      name = "testActorFirstBytes",
    )

    val src = DefaultPathBuilder.createTempFile()
    src.write("hello")

    val readCommand = DefaultIoContentAsStringCommand(src, IoReadOptions(Option(2), failOnOverflow = false))

    testActor ! readCommand
    expectMsgPF(5 seconds) {
      case response: IoSuccess[_] =>
        response.command.isInstanceOf[IoContentAsStringCommand] shouldBe true
        response.result.asInstanceOf[String] shouldBe "he"
      case response: IoFailure[_] => fail("Expected an IoSuccess", response.failure)
    }

    src.delete()
  }

  it should "read the file if it's under the byte limit" in {
    val testActor = TestActorRef(
      factory = new IoActor(1, 10, 10, None, TestProbe("serviceRegistryActorByteLimit").ref, "cromwell test"),
      name = "testActorByteLimit",
    )

    val src = DefaultPathBuilder.createTempFile()
    src.write("hello")

    val readCommand = DefaultIoContentAsStringCommand(src, IoReadOptions(Option(6), failOnOverflow = true))

    testActor ! readCommand
    expectMsgPF(5 seconds) {
      case response: IoSuccess[_] =>
        response.command.isInstanceOf[IoContentAsStringCommand] shouldBe true
        response.result.asInstanceOf[String] shouldBe "hello"
      case response: IoFailure[_] => fail("Expected an IoSuccess", response.failure)
    }

    src.delete()
  }

  it should "fail if the file is larger than the read limit" in {
    val testActor = TestActorRef(
      factory = new IoActor(1, 10, 10, None, TestProbe("serviceRegistryActorReadLimit").ref, "cromwell test"),
      name = "testActorReadLimit",
    )

    val src = DefaultPathBuilder.createTempFile()
    src.write("hello")

    val readCommand = DefaultIoContentAsStringCommand(src, IoReadOptions(Option(2), failOnOverflow = true))

    testActor ! readCommand
    expectMsgPF(5 seconds) {
      case _: IoSuccess[_] => fail("Command should have failed because the read limit was < file size and failOnOverflow was true")
      case response: IoFailure[_] => response.failure.getMessage shouldBe s"[Attempted 1 time(s)] - IOException: Could not read from ${src.pathAsString}: File ${src.pathAsString} is larger than requested maximum of 2 Bytes."
    }

    src.delete()
  }

  it should "return a file size" in {
    val testActor = TestActorRef(
      factory = new IoActor(1, 10, 10, None, TestProbe("serviceRegistryActorSize").ref, "cromwell test"),
      name = "testActorSize",
    )

    val src = DefaultPathBuilder.createTempFile()
    src.write("hello")

    val sizeCommand = DefaultIoSizeCommand(src)

    testActor ! sizeCommand
    expectMsgPF(5 seconds) {
      case response: IoSuccess[_] =>
        response.command.isInstanceOf[IoSizeCommand] shouldBe true
        response.result.asInstanceOf[Long] shouldBe 5
      case response: IoFailure[_] => fail("Expected an IoSuccess", response.failure)
    }

    src.delete()
  }

  it should "return a file md5 hash (local)" in {
    val testActor = TestActorRef(
      factory = new IoActor(1, 10, 10, None, TestProbe("serviceRegistryActorHash").ref, "cromwell test"),
      name = "testActorHash",
    )

    val src = DefaultPathBuilder.createTempFile()
    src.write("hello")

    val hashCommand = DefaultIoHashCommand(src)

    testActor ! hashCommand
    expectMsgPF(5 seconds) {
      case response: IoSuccess[_] =>
        response.command.isInstanceOf[IoHashCommand] shouldBe true
        response.result.asInstanceOf[String] shouldBe "5d41402abc4b2a76b9719d911017c592"
      case response: IoFailure[_] => fail("Expected an IoSuccess", response.failure)
    }

    src.delete()
  }

  it should "touch a file (local)" in {
    val testActor = TestActorRef(
      factory = new IoActor(1, 10, 10, None, TestProbe("serviceRegistryActorTouch").ref, "cromwell test"),
      name = "testActorTouch",
    )

    val src = DefaultPathBuilder.createTempFile()
    src.write("hello")

    val touchCommand = DefaultIoTouchCommand(src)

    testActor ! touchCommand
    expectMsgPF(5 seconds) {
      case _: IoSuccess[_] =>
      case response: IoFailure[_] => fail("Expected an IoSuccess", response.failure)
    }

    src.delete()
  }

  it should "have correct retryable exceptions" in {
    val retryables = List(
      new StorageException(500, "message"),
      new StorageException(502, "message"),
      new StorageException(503, "message"),
      new StorageException(504, "message"),
      new StorageException(408, "message"),
      new StorageException(429, "message"),
      BatchFailedException(new Exception),
      new SocketException(),
      new SocketTimeoutException(),
      new IOException("text Error getting access token for service account some other text"),

      new IOException("Could not read from gs://fc-secure-<snip>/JointGenotyping/<snip>/call-HardFilterAndMakeSitesOnlyVcf/shard-4688/rc: 500 Internal Server Error\nBackend Error"),
      new IOException("Could not read from gs://fc-secure-<snip>/JointGenotyping/<snip>/call-HardFilterAndMakeSitesOnlyVcf/shard-4688/rc: 500 Internal Server Error Backend Error"),

      new IOException("Could not read from gs://broad-epi-cromwell/workflows/ChipSeq/ce6a5671-baf6-4734-a32b-abf3d9138e9b/call-epitope_classifier/memory_retry_rc: 503 Service Unavailable\nBackend Error"),
      new IOException("Could not read from gs://fc-secure-<snip>/JointGenotyping/<snip>/call-HardFilterAndMakeSitesOnlyVcf/shard-4688/rc: 503 Service Unavailable Backend Error"),

      new IOException("Could not read from gs://mccarroll-mocha/cromwell/cromwell-executions/mocha/86d47e9a-5745-4ec0-b4eb-0164f073e5f4/call-idat2gtc/shard-73/rc: 504 Gateway Timeout\nGET https://storage.googleapis.com/download/storage/v1/b/mccarroll-mocha/o/cromwell%2Fcromwell-executions%2Fmocha%2F86d47e9a-5745-4ec0-b4eb-0164f073e5f4%2Fcall-idat2gtc%2Fshard-73%2Frc?alt=media"),

      // Prove that `isRetryable` successfully recurses to unwrap the lowest-level Throwable
      new IOException(new Throwable("Could not read from gs://fc-secure-<snip>/JointGenotyping/<snip>/call-HardFilterAndMakeSitesOnlyVcf/shard-4688/rc: 500 Internal Server Error Backend Error")),
      new IOException(new Throwable("Could not read from gs://fc-secure-<snip>/JointGenotyping/<snip>/call-HardFilterAndMakeSitesOnlyVcf/shard-4688/rc: 503 Service Unavailable Backend Error")),

      new IOException("Some other text. Could not read from gs://fc-secure-<snip>/JointGenotyping/<snip>/call-HardFilterAndMakeSitesOnlyVcf/shard-4688/rc: 503 Service Unavailable"),
      new IOException("Some other text. Could not read from gs://fc-secure-<snip>/JointGenotyping/<snip>/call-HardFilterAndMakeSitesOnlyVcf/shard-4688/rc: 504 Gateway Timeout"),
    )

    retryables foreach { RetryableRequestSupport.isRetryable(_) shouldBe true }
  }

  it should "have correct non-retryable exceptions" in {
    val nonRetryables = List(
      new IOException("message: 500 Internal Server Error Backend Error"),
      new IOException("404 File Not Found"),
      new IOException("502 HTTP Status Code"),
      new Exception("502 HTTP Status Code"),
      new Exception("5xx HTTP Status Code"),
      new IOException("Could not read from gs://fc-secure-<snip>/JointGenotyping/<snip>/call-HardFilterAndMakeSitesOnlyVcf/shard-500/rc: 404 File Not Found")
    )

    nonRetryables foreach { RetryableRequestSupport.isRetryable(_) shouldBe false }
  }

  it should "not crash when certain exception members are `null`" in {

    // Javadoc for `com.google.cloud.storage.StorageException` says `message`, `cause` may be `null`
    val nullCause = new StorageException(3, "blah", "no reason", null)
    val nullMessage = new StorageException(4, null)

    RetryableRequestSupport.isRetryable(nullCause) shouldBe false
    RetryableRequestSupport.isRetryable(nullMessage) shouldBe false

  }
}
