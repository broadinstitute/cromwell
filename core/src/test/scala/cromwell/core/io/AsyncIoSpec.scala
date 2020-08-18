package cromwell.core.io

import java.nio.file.{FileAlreadyExistsException, NoSuchFileException}
import java.util.UUID

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.testkit.TestActorRef
import cromwell.core.TestKitSuite
import cromwell.core.path.DefaultPathBuilder
import org.scalatest.flatspec.AsyncFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatestplus.mockito.MockitoSugar

class AsyncIoSpec extends TestKitSuite with AsyncFlatSpecLike with Matchers with MockitoSugar {

  behavior of "AsyncIoSpec"
  
  implicit val ioCommandBuilder = DefaultIoCommandBuilder
  
  it should "write asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()
    
    testActor.underlyingActor.asyncIo.writeAsync(testPath, "hello", Seq.empty) map { _ =>
      assert(testPath.contentAsString == "hello")
    }
  }

  it should "read asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("hello")

    testActor.underlyingActor.asyncIo.contentAsStringAsync(testPath, None, failOnOverflow = false) map { result =>
      assert(result == "hello")
    }
  }

  it should "get size asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("hello")

    testActor.underlyingActor.asyncIo.sizeAsync(testPath) map { size =>
      assert(size == 5)
    }
  }

  it should "get hash asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("hello")

    testActor.underlyingActor.asyncIo.hashAsync(testPath) map { hash =>
      assert(hash == "5D41402ABC4B2A76B9719D911017C592")
    }
  }

  it should "copy asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()
    val testCopyPath = testPath.sibling(UUID.randomUUID().toString)

    testActor.underlyingActor.asyncIo.copyAsync(testPath, testCopyPath) map { hash =>
      assert(testCopyPath.exists)
    }

    testPath.write("new text")
    
    // Honor overwrite true
    testActor.underlyingActor.asyncIo.copyAsync(testPath, testCopyPath, overwrite = true) map { hash =>
      assert(testCopyPath.exists)
      assert(testCopyPath.contentAsString == "new text")
    }

    // Honor overwrite false
    recoverToSucceededIf[FileAlreadyExistsException] { testActor.underlyingActor.asyncIo.copyAsync(testPath, testCopyPath, overwrite = false) }
  }

  it should "delete asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()

    testActor.underlyingActor.asyncIo.deleteAsync(testPath) map { _ =>
      assert(!testPath.exists)
    }

    // Honor swallow exception true
    testActor.underlyingActor.asyncIo.deleteAsync(testPath, swallowIoExceptions = true) map { _ =>
      assert(!testPath.exists)
    }

    // Honor swallow exception false
    recoverToSucceededIf[NoSuchFileException] { testActor.underlyingActor.asyncIo.deleteAsync(testPath, swallowIoExceptions = false) }
  }

  private class AsyncIoTestActor(override val ioActor: ActorRef) extends Actor with ActorLogging with AsyncIoActorClient {

    override lazy val ioCommandBuilder = DefaultIoCommandBuilder
    
    override def receive: Receive = {
      case _ =>
    }
  }

}
