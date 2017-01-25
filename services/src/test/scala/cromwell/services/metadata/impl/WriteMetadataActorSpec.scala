package cromwell.services.metadata.impl

import java.time.OffsetDateTime

import akka.testkit.TestFSMRef
import cats.data.NonEmptyVector
import cromwell.core.WorkflowId
import cromwell.services.ServicesSpec
import cromwell.services.metadata.MetadataService.PutMetadataAction
import cromwell.services.metadata.{MetadataEvent, MetadataKey, MetadataValue}
import cromwell.services.metadata.impl.WriteMetadataActor.{HasEvents, NoEvents, WaitingToWrite, WritingToDb}

import scala.concurrent.duration._

class WriteMetadataActorSpec extends ServicesSpec("Metadata") {
  import WriteMetadataActorSpec.Action

  val actor = TestFSMRef(WriteMetadataActor(10, 1.days)) // A WMA that won't (hopefully!) perform a time based flush during this test

  "WriteMetadataActor" should {
    "start with no events and waiting to write" in {
      assert(actor.stateName == WaitingToWrite)
      assert(actor.stateData == NoEvents)
    }

    "Have one event and be waiting after one event is sent" in {
      actor ! Action
      assert(actor.stateName == WaitingToWrite)
      assert(actor.stateData == HasEvents(NonEmptyVector.fromVectorUnsafe(Action.events.toVector)))
    }

    "Have one event after batch size + 1 is reached" in {
      1 to 10 foreach { _ => actor ! Action }
      assert(actor.stateName == WritingToDb)
      actor ! Action
      assert(actor.stateData == HasEvents(NonEmptyVector.fromVectorUnsafe(Action.events.toVector)))
    }
  }
}

object WriteMetadataActorSpec {
  val Event = MetadataEvent(MetadataKey(WorkflowId.randomId(), None, "key"), Option(MetadataValue("value")), OffsetDateTime.now)
  val Action = PutMetadataAction(Event)
}
