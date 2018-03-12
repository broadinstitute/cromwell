package cromwell.engine.io

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
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._
import scala.language.postfixOps

class IoActorSpec extends TestKitSuite with FlatSpecLike with Matchers with ImplicitSender {
  behavior of "IoActor"
  
  implicit val actorSystem = system
  implicit val ec: ExecutionContext = system.dispatcher
  implicit val materializer = ActorMaterializer()
  
  override def afterAll() = {
    materializer.shutdown()
    super.afterAll()
  }
  
  it should "copy a file" in {
    val testActor = TestActorRef(new IoActor(1, 10, 10, None, TestProbe().ref))
    
    val src = DefaultPathBuilder.createTempFile()
    val dst: Path = src.parent.resolve(src.name + "-dst")
    
    val copyCommand = DefaultIoCopyCommand(src, dst, overwrite = true)
    
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
    val testActor = TestActorRef(new IoActor(1, 10, 10, None, TestProbe().ref))

    val src = DefaultPathBuilder.createTempFile()

    val writeCommand = DefaultIoWriteCommand(src, "hello", OpenOptions.default)

    testActor ! writeCommand
    expectMsgPF(5 seconds) {
      case response: IoSuccess[_] => response.command.isInstanceOf[IoWriteCommand] shouldBe true
      case response: IoFailure[_] => fail("Expected an IoSuccess", response.failure)
    }

    src.contentAsString shouldBe "hello"
    src.delete()
  }

  it should "delete a file" in {
    val testActor = TestActorRef(new IoActor(1, 10, 10, None, TestProbe().ref))

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
    val testActor = TestActorRef(new IoActor(1, 10, 10, None, TestProbe().ref))

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
    val testActor = TestActorRef(new IoActor(1, 10, 10, None, TestProbe().ref))

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

  it should "fail if the file is larger than the read limit" in {
    val testActor = TestActorRef(new IoActor(1, 10, 10, None, TestProbe().ref))

    val src = DefaultPathBuilder.createTempFile()
    src.write("hello")

    val readCommand = DefaultIoContentAsStringCommand(src, IoReadOptions(Option(2), failOnOverflow = true))

    testActor ! readCommand
    expectMsgPF(5 seconds) {
      case _: IoSuccess[_] => fail("Command should have failed because the read limit was < file size and failOnOverflow was true")
      case response: IoFailure[_] => response.failure.getMessage shouldBe s"java.io.IOException: File ${src.pathAsString} is larger than 2 Bytes"
    }

    src.delete()
  }

  it should "return a file size" in {
    val testActor = TestActorRef(new IoActor(1, 10, 10, None, TestProbe().ref))

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
    val testActor = TestActorRef(new IoActor(1, 10, 10, None, TestProbe().ref))

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
    val testActor = TestActorRef(new IoActor(1, 10, 10, None, TestProbe().ref))

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
      new SocketTimeoutException()
    )

    retryables foreach { IoActor.isRetryable(_) shouldBe true }
    retryables foreach { IoActor.isFatal(_) shouldBe false }
  }
}
