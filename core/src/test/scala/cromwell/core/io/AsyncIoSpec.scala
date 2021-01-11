package cromwell.core.io

import java.nio.file.NoSuchFileException
import java.util.UUID

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.testkit.TestActorRef
import cromwell.core.TestKitSuite
import cromwell.core.path.{DefaultPathBuilder, Path}
import org.scalatest.flatspec.AsyncFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatestplus.mockito.MockitoSugar

import scala.util.{Failure, Try}
import scala.util.control.NoStackTrace

class AsyncIoSpec extends TestKitSuite with AsyncFlatSpecLike with Matchers with MockitoSugar {

  behavior of "AsyncIoSpec"
  
  implicit val ioCommandBuilder: DefaultIoCommandBuilder.type = DefaultIoCommandBuilder
  
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

    testActor.underlyingActor.asyncIo.copyAsync(testPath, testCopyPath) map { _ =>
      assert(testCopyPath.exists)
    }

    testPath.write("new text")
    
    testActor.underlyingActor.asyncIo.copyAsync(testPath, testCopyPath) map { _ =>
      assert(testCopyPath.exists)
      assert(testCopyPath.contentAsString == "new text")
    }
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
    //noinspection RedundantDefaultArgument
    recoverToSucceededIf[NoSuchFileException] {
      testActor.underlyingActor.asyncIo.deleteAsync(testPath, swallowIoExceptions = false)
    }
  }

  it should "handle command creation errors asynchronously" in {
    val partialIoCommandBuilder = new PartialIoCommandBuilder {
      override def existsCommand: PartialFunction[Path, Try[IoExistsCommand]] = {
        case _ => Failure(new Exception("everything's fine, I am an expected exists fail") with NoStackTrace)
      }
    }
    val testActor =
      TestActorRef(new AsyncIoTestActor(simpleIoActor, new IoCommandBuilder(List(partialIoCommandBuilder))))

    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("hello")

    testActor.underlyingActor.asyncIo.existsAsync(testPath).failed map { throwable =>
      assert(throwable.getMessage == "everything's fine, I am an expected exists fail")
    }
  }

  private class AsyncIoTestActor(override val ioActor: ActorRef,
                                 override val ioCommandBuilder: IoCommandBuilder = DefaultIoCommandBuilder
                                ) extends Actor with ActorLogging with AsyncIoActorClient {

    override def receive: Receive = {
      case _ =>
    }
  }

}
