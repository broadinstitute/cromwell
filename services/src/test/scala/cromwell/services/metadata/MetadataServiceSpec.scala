package cromwell.services.metadata

import cromwell.core.WorkflowId
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

}
