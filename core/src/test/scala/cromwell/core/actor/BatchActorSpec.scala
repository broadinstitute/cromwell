package cromwell.core.actor

import akka.actor.ActorRef
import akka.routing.Listen
import akka.testkit.{TestFSMRef, TestProbe}
import cats.data.NonEmptyVector
import cromwell.core.TestKitSuite
import cromwell.core.actor.BatchActor._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.scalatest.concurrent.Eventually
import org.scalatest.{FlatSpecLike, Matchers}

import scala.concurrent.duration._
import scala.concurrent.{Future, Promise}

class BatchActorSpec extends TestKitSuite with FlatSpecLike with Matchers with Eventually {

  behavior of "BatchingDbWriter"
  
  override val patienceConfig = PatienceConfig(timeout = scaled(5.seconds), interval = scaled(1.second))
  implicit val patience = patienceConfig

  it should "start with WaitingToProcess" in {
    val batch = TestFSMRef(new BatchActorTest)
    batch.stateName shouldBe WaitingToProcess
  }

  it should "enqueue commands" in {
    val batch = TestFSMRef(new BatchActorTest)
    batch ! "hello"
    batch.stateData.weight shouldBe 5
    batch.stateData.innerQueue.size shouldBe 1
  }

  it should "immediately process the head if and only if the next messages puts the queue weight over batch size" in {
    val batch = TestFSMRef(new BatchActorTest)
    batch ! "hello"
    batch ! "hola"
    batch.underlyingActor.processed shouldBe Vector.empty
    // With this message the weight goes over 10
    batch ! "bonjour"
    
    batch.underlyingActor.processed shouldBe Vector("hello", "hola")

    batch.stateData.weight shouldBe 7
    batch.stateData.innerQueue.size shouldBe 1
    batch.stateData.innerQueue.head shouldBe "bonjour"
  }

  it should "process the head of the queue when receiving a ScheduledProcessAction and broadcast its queue weight" in {
    val batch = TestFSMRef(new BatchActorTest)
    val listener = TestProbe()
    batch ! Listen(listener.ref)
    batch ! "hello"
    batch ! "hola"

    batch.stateData.weight shouldBe 9
    batch.stateData.innerQueue.size shouldBe 2

    batch ! ScheduledProcessAction
    listener.expectMsgPF(5.second) {
      case qw: QueueWeight => qw.weight shouldBe 2
      case _ =>
    }
    batch.underlyingActor.processed shouldBe Vector("hello", "hola")
    batch.stateData.weight shouldBe 0
    batch.stateData.innerQueue.size shouldBe 0
  }

  it should "keep enqueueing events in Processing state" in {
    val batch = TestFSMRef(new BatchActorTest)
    batch ! "hello"
    batch.setState(Processing)
    batch ! "hola"

    batch.stateData.weight shouldBe 9
    batch.stateData.innerQueue.size shouldBe 2
  }

  it should "stay where it is when receiving ScheduledProcessAction in Processing state and broadcast its queue weight" in {
    val batch = TestFSMRef(new BatchActorTest)
    val listener = TestProbe()
    batch ! Listen(listener.ref)
    batch.setState(Processing)
    batch ! "hola"

    batch.stateData.weight shouldBe 4
    batch.stateData.innerQueue.size shouldBe 1
    batch ! ScheduledProcessAction
    listener.expectMsgPF(5.second) {
      case qw: QueueWeight => qw.weight shouldBe 1
      case _ =>
    }
  }

  it should "go back to WaitingToProcess after processing if not over batch size" in {
    // Add some processing time so we have time to send more events while a batch is being processed
    val batch = TestFSMRef(new BatchActorTest(2.seconds))
    batch ! "hola"
    batch ! "hello"
    batch ! "bonjour"

    // at that point we've processed the first batch and are waiting for the future to complete
    batch.underlyingActor.processed shouldBe Vector("hola", "hello")

    eventually {
      batch.stateName shouldBe WaitingToProcess
    }
  }

  it should "process again when previous processing finished and we're still over batch size" in {
    // Add some processing time so we have time to send more events while a batch is being processed
    val batch = TestFSMRef(new BatchActorTest(2.seconds))
    batch ! "hola"
    batch ! "hello"
    batch ! "bonjour"

    // at that point we've processed the first batch and are waiting for the future to complete
    batch.underlyingActor.processed shouldBe Vector("hola", "hello")

    // send more messages to stay over the batch size when the future completes
    batch ! "aurevoir"
    batch ! "goodbye"

    eventually {
      batch.underlyingActor.processed shouldBe Vector("hola", "hello", "bonjour", "aurevoir")
    }
  }

  it should "shutdown when in WaitingToProcess if queue is empty" in {
    val batch = TestFSMRef(new BatchActorTest)
    val watcher = TestProbe()
    watcher.watch(batch)
    batch ! ShutdownCommand
    
    eventually {
      watcher.expectTerminated(batch)
    }
  }

  it should "flush and shutdown when in WaitingToProcess if queue is not empty" in {
    // Add some processing time so we have time to check elements have been processed before the actor dies
    val batch = TestFSMRef(new BatchActorTest(2.seconds))
    val watcher = TestProbe()
    watcher.watch(batch)
    batch ! "hola"
    batch ! "hello"
    batch ! ShutdownCommand

    eventually {
      batch.underlyingActor.processed == Vector("hola", "hello")
      watcher.expectTerminated(batch)
    }
  }

  it should "flush everything and stop if received a shutdown command when Processing" in {
    // Add some processing time so we have time to check elements have been processed before the actor dies
    val batch = TestFSMRef(new BatchActorTest(2.seconds))
    batch.setState(Processing)
    val watcher = TestProbe()
    watcher.watch(batch)
    batch ! "hola"
    batch ! "hello"
    
    batch ! ShutdownCommand
    
    batch ! ProcessingComplete

    eventually {
      batch.underlyingActor.processed == Vector("hola", "hello")
      watcher.expectTerminated(batch)
    }
  }

  it should "keep going even if processing fail" in {
    // Add some processing time so we have time to send more events while a batch is being processed
    val batch = TestFSMRef(new BatchActorTest(Duration.Zero, true))
    batch ! "bonjour"
    batch ! "hello"

    batch.underlyingActor.processed shouldBe Vector.empty

    eventually {
      // Even if writing bonjour fails we go back to WaitingToProcess and "hello" is still there to be processed later
      batch.stateName shouldBe WaitingToProcess
      batch.stateData.weight shouldBe 5
    }
  }

  class BatchActorTest(processingTime: FiniteDuration = Duration.Zero, fail: Boolean = false) extends BatchActor[String](10.hours, 10) {
    var processed: Vector[String] = Vector.empty
    override def commandToData(snd: ActorRef) = {
      case command: String => command
    }
    override protected def weightFunction(command: String) = command.length
    override protected def process(data: NonEmptyVector[String]) = {
      if (processingTime != Duration.Zero) {
        processed = processed ++ data.toVector
       val promise = Promise[Int]
        system.scheduler.scheduleOnce(processingTime) { promise.success(data.map(weightFunction).toVector.sum) }
        promise.future
      } else if (!fail) {
        processed = processed ++ data.toVector
        Future.successful(data.map(weightFunction).toVector.sum)
      }
      else Future.failed(new Exception("Oh nose ! (This is a test failure and is expected !)"))
    }
  }

}
