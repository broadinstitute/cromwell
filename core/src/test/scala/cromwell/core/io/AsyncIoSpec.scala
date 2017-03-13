package cromwell.core.io

import java.nio.file.{FileAlreadyExistsException, NoSuchFileException}
import java.util.UUID

import akka.actor.{Actor, ActorLogging, ActorRef}
import akka.testkit.TestActorRef
import cromwell.core.path.DefaultPathBuilder
import cromwell.core.{SimpleIoActor, TestKitSuite}
import org.scalatest.mockito.MockitoSugar
import org.scalatest.{AsyncFlatSpecLike, Matchers}

class AsyncIoSpec extends TestKitSuite with AsyncFlatSpecLike with Matchers with MockitoSugar {

  behavior of "AsyncIoSpec"
  
  val simpleIoActor = system.actorOf(SimpleIoActor.props)
  
  override def afterAll() = {
    system stop simpleIoActor
    super.afterAll()
  } 
  
  it should "write asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()
    
    testActor.underlyingActor.writeAsync(testPath, "hello", Seq.empty) map { _ =>
      assert(testPath.contentAsString == "hello")
    }
  }

  it should "read asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("hello")

    testActor.underlyingActor.contentAsStringAsync(testPath) map { result =>
      assert(result == "hello")
    }
  }

  it should "get size asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("hello")

    testActor.underlyingActor.sizeAsync(testPath) map { size =>
      assert(size == 5)
    }
  }

  it should "get hash asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()
    testPath.write("hello")

    testActor.underlyingActor.hashAsync(testPath) map { hash =>
      assert(hash == "5D41402ABC4B2A76B9719D911017C592")
    }
  }

  it should "copy asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()
    val testCopyPath = testPath.sibling(UUID.randomUUID().toString)

    testActor.underlyingActor.copyAsync(testPath, testCopyPath) map { hash =>
      assert(testCopyPath.exists)
    }

    testPath.write("new text")
    
    // Honor overwrite true
    testActor.underlyingActor.copyAsync(testPath, testCopyPath, overwrite = true) map { hash =>
      assert(testCopyPath.exists)
      assert(testCopyPath.contentAsString == "new text")
    }

    // Honor overwrite false
    recoverToSucceededIf[FileAlreadyExistsException] { testActor.underlyingActor.copyAsync(testPath, testCopyPath, overwrite = false) }
  }

  it should "delete asynchronously" in {
    val testActor = TestActorRef(new AsyncIoTestActor(simpleIoActor))

    val testPath = DefaultPathBuilder.createTempFile()

    testActor.underlyingActor.deleteAsync(testPath) map { _ =>
      assert(!testPath.exists)
    }

    // Honor swallow exception true
    testActor.underlyingActor.deleteAsync(testPath, swallowIoExceptions = true) map { _ =>
      assert(!testPath.exists)
    }

    // Honor swallow exception false
    recoverToSucceededIf[NoSuchFileException] { testActor.underlyingActor.deleteAsync(testPath, swallowIoExceptions = false) }
  }

  private class AsyncIoTestActor(override val ioActor: ActorRef) extends Actor with ActorLogging with AsyncIo with DefaultIoCommandBuilder {

    context.become(ioReceive orElse receive)

    override def receive: Receive = {
      case _ =>
    }
  }

}
