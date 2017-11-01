package cromwell.services.metadata.impl

import java.time.OffsetDateTime
import java.util.UUID

import akka.testkit.TestFSMRef
import cats.data.NonEmptyVector
import cromwell.core.WorkflowId
import cromwell.core.actor.BatchingDbWriter
import cromwell.core.actor.BatchingDbWriter._
import cromwell.services.ServicesSpec
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import org.scalatest.BeforeAndAfter
import org.scalatest.concurrent.Eventually

import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.Success

class WriteMetadataActorSpec extends ServicesSpec("Metadata") with Eventually with BeforeAndAfter {
  import WriteMetadataActorSpec.Action

  var actor: TestFSMRef[BatchingDbWriterState, BatchingDbWriter.BatchingDbWriterData, DelayingWriteMetadataActor] = _

  before {
    actor = TestFSMRef(new DelayingWriteMetadataActor(), "DelayingWriteMetadataActor-" + UUID.randomUUID())
  }

  "WriteMetadataActor" should {
    "start with no events and waiting to write" in {
      actor.stateName shouldBe WaitingToWrite
      actor.stateData shouldBe NoData
    }

    "Have one event and be waiting after one event is sent" in {
      actor ! Action
      eventually {
        actor.stateName shouldBe WaitingToWrite
        actor.stateData shouldBe HasData(NonEmptyVector.fromVectorUnsafe(Action.events.toVector))
      }
    }

    "Have one event after batch size + 1 is reached" in {
      1 to WriteMetadataActorSpec.BatchRate foreach { _ => actor ! Action }
      actor.stateName shouldBe WaitingToWrite

      eventually {
        actor.stateData match {
          case HasData(e) => e.toVector.size shouldBe WriteMetadataActorSpec.BatchRate
          case _ => fail("Expecting the actor to have events queued up")
        }
      }
      actor ! Action
      eventually {
        actor.stateName shouldBe WritingToDb
        actor.underlyingActor.writeToDbInProgress shouldBe true
        actor.stateData shouldBe NoData
      }
      actor ! Action
      eventually {
        actor.stateName shouldBe WritingToDb
        actor.stateData shouldBe HasData(NonEmptyVector.fromVectorUnsafe(Action.events.toVector))
      }
      actor.underlyingActor.completeWritePromise()
      eventually {
        actor.stateName shouldBe WaitingToWrite
        actor.stateData shouldBe HasData(NonEmptyVector.fromVectorUnsafe(Action.events.toVector))
      }
    }
  }
}

object WriteMetadataActorSpec {
  val Event = MetadataEvent(MetadataKey(WorkflowId.randomId(), None, "key"), Option(MetadataValue("value")), OffsetDateTime.now)
  val Action = PutMetadataAction(Event)

  val BatchRate: Int = 10
  val FunctionallyForever: FiniteDuration = 100.days
}

// A WMA that won't (hopefully!) perform a time based flush during this test
final class DelayingWriteMetadataActor extends WriteMetadataActor(WriteMetadataActorSpec.BatchRate, WriteMetadataActorSpec.FunctionallyForever) {

  var writeToDbInProgress: Boolean = false
  var writeToDbCompletionPromise: Option[Promise[Unit]] = None

  override def addMetadataEvents(metadataEvents: Iterable[MetadataEvent])(implicit ec: ExecutionContext): Future[Unit] = {
    writeToDbCompletionPromise = Option(Promise[Unit]())
    writeToDbInProgress = true
    writeToDbCompletionPromise.get.future
  }

  def completeWritePromise(): Unit = {
    writeToDbCompletionPromise match {
      case Some(promise) =>
        promise.complete(Success(()))
        writeToDbInProgress = false
        writeToDbCompletionPromise = None
      case None => throw new Exception("BAD TEST! Cannot complete the actor's write future if the actor hasn't requested it yet!")
    }
  }
}
