package cromwell.services.metadata

import java.util.UUID

import cromwell.core.WorkflowId
import lenthall.exception.AggregatedException
import org.scalactic.Equality
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wdl4s.types.{WdlArrayType, WdlMapType, WdlStringType}
import wdl4s.values._

class MetadataServiceSpec extends FlatSpec with Matchers with TableDrivenPropertyChecks {

  behavior of "MetadataServiceSpec"

  // Discard timestamp when determining equality here
  implicit val metadataEventEquality = new Equality[MetadataEvent] {
    override def areEqual(a: MetadataEvent, b: Any): Boolean = b match {
      case bEvent: MetadataEvent =>
        a.key == bEvent.key && a.value == bEvent.value
      case _ => false
    }
  }

  it should "convert a WdlArray to MetadataEvents" in {
    import MetadataService._
    val workflowId = WorkflowId.randomId()
    val wdlArray = WdlArray(WdlArrayType(WdlStringType), Seq(WdlString("Hello"), WdlString("world!")))
    val emptyWdlArray = WdlArray(WdlArrayType(WdlStringType), Seq.empty)

    wdlValueToMetadataEvents(MetadataKey(workflowId, None, "root"), wdlArray).toList should contain theSameElementsInOrderAs List(
      MetadataEvent(MetadataKey(workflowId, None, "root[0]"), MetadataValue("Hello")),
      MetadataEvent(MetadataKey(workflowId, None, "root[1]"), MetadataValue("world!"))
    )

    wdlValueToMetadataEvents(MetadataKey(workflowId, None, "root"), emptyWdlArray).toList should contain theSameElementsAs List(
      MetadataEvent.empty(MetadataKey(workflowId, None, "root[]"))
    )
  }

  it should "convert a WdlMap to MetadataEvents" in {
    import MetadataService._
    val workflowId = WorkflowId.randomId()
    val wdlArray = WdlMap(WdlMapType(WdlStringType, WdlStringType), Map(
      WdlString("Hello") -> WdlString("world!"),
      WdlString("Goodbye") -> WdlString("world!")
    ))
    val emptyWdlMap = WdlMap(WdlMapType(WdlStringType, WdlStringType), Map.empty)

    wdlValueToMetadataEvents(MetadataKey(workflowId, None, "root"), wdlArray).toList should contain theSameElementsInOrderAs List(
      MetadataEvent(MetadataKey(workflowId, None, "root:Hello"), MetadataValue("world!")),
      MetadataEvent(MetadataKey(workflowId, None, "root:Goodbye"), MetadataValue("world!"))
    )

    wdlValueToMetadataEvents(MetadataKey(workflowId, None, "root"), emptyWdlMap).toList should contain theSameElementsAs List(
      MetadataEvent.empty(MetadataKey(workflowId, None, "root"))
    )
  }

  it should "convert a primitive WdlValue to MetadataEvents" in {
    import MetadataService._
    val workflowId = WorkflowId.randomId()

    val values = Table(
      ("wdlValue", "metadataValue"),
      (WdlString("hi"), MetadataValue("hi", MetadataString)),
      (WdlInteger(1), MetadataValue("1", MetadataInt)),
      (WdlFloat(1F), MetadataValue("1.0", MetadataNumber)),
      (WdlBoolean(true), MetadataValue("true", MetadataBoolean))
    )

    forAll(values) { (wdlValue, metadataValue) =>
      wdlValueToMetadataEvents(MetadataKey(workflowId, None, "root"), wdlValue).toList should contain theSameElementsAs List(
        MetadataEvent(MetadataKey(workflowId, None, "root"), metadataValue)
      )
    }
  }

  // For the metadata event tests!
  val eventKeyRegex = "path\\:to\\:failures\\[([0-9]*)\\]\\:message".r

  it should "convert an exception into a failure event" in {
    import MetadataService.throwableToMetadataEvents

    val workflowId = WorkflowId(UUID.randomUUID())
    val mdkey = MetadataKey(workflowId, None, "path:to:failures")

    val tMsg = "The Oscars suck!"
    val t = new RuntimeException(tMsg)

    val events = throwableToMetadataEvents(mdkey, t)
    events.size should be(1)
    val _ = validateEvent(events.head, workflowId, tMsg)
  }

  it should "convert nested exceptions into a sequence of failure events" in {
    import MetadataService.throwableToMetadataEvents

    val workflowId = WorkflowId(UUID.randomUUID())
    val mdkey = MetadataKey(workflowId, None, "path:to:failures")

    val innerCauseMsg = "Envelope malfunctions"
    val innerCause = new RuntimeException(innerCauseMsg)

    val causeMsg = "Wrong recipients"
    val cause = new RuntimeException(causeMsg, innerCause)

    val tMsg = "The Oscars suck!"
    val t = new RuntimeException(tMsg, cause)

    val events = throwableToMetadataEvents(mdkey, t)
    events.size should be(3)
    val errorIds = List(
      validateEvent(events.head, workflowId, tMsg),
      validateEvent(events(1), workflowId, causeMsg),
      validateEvent(events(2), workflowId, innerCauseMsg)
    )
    // Make sure all IDs are unique:
    errorIds.toSet.size should be(3)
  }

  it should "convert aggregated exceptions into a sequence of failure events" in {
    import MetadataService.throwableToMetadataEvents

    val workflowId = WorkflowId(UUID.randomUUID())
    val mdkey = MetadataKey(workflowId, None, "path:to:failures")

    val innerCauseMsg = "Envelope malfunctions"
    val innerCause = new RuntimeException(innerCauseMsg)

    val causeMsg1 = "Wrong recipients"
    val cause1 = new RuntimeException(causeMsg1, innerCause)
    val causeMsg2 = "Self congratulation"
    val cause2 = new RuntimeException(causeMsg2)
    val causeMsg3 = "The Globes are better anyway"
    val cause3 = new RuntimeException(causeMsg3)

    val causeContext = "Compound Entertainment Failure"
    val cause = new AggregatedException(causeContext, List(cause1, cause2, cause3))

    val tMsg = "The Oscars suck!"
    val t = new RuntimeException(tMsg, cause)

    val events = throwableToMetadataEvents(mdkey, t)
    events.size should be(6)
    val errorIds = List(
      validateEvent(events.head, workflowId, tMsg),
      validateEvent(events(1), workflowId, causeContext),
      validateEvent(events(2), workflowId, causeMsg1),
      validateEvent(events(3), workflowId, innerCauseMsg),
      validateEvent(events(4), workflowId, causeMsg2),
      validateEvent(events(5), workflowId, causeMsg3)
    )
    // Make sure all IDs are unique:
    errorIds.toSet.size should be(6)
  }

  def validateEvent(event: MetadataEvent, workflowId: WorkflowId, message: String) = event match {
    case MetadataEvent(k, Some(MetadataValue(v, _)), _) =>
      k.workflowId should be(workflowId)
      val messageIndex = k.key match {
        case eventKeyRegex(x) => x
        case _ => fail("Unexpected failure key format: " + k.key)
      }
      v should be(message)

      messageIndex
  }
}
