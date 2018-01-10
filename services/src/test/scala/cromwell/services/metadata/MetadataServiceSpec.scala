package cromwell.services.metadata

import java.util.UUID

import cromwell.core.WorkflowId
import common.exception.AggregatedException
import org.scalactic.Equality
import org.scalatest.prop.TableDrivenPropertyChecks
import org.scalatest.{FlatSpec, Matchers}
import wom.types._
import wom.values._

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

  it should "convert a WomArray to MetadataEvents" in {
    import MetadataService._
    val workflowId = WorkflowId.randomId()
    val womArray = WomArray(WomArrayType(WomStringType), Seq(WomString("Hello"), WomString("world!")))
    val emptyWomArray = WomArray(WomArrayType(WomStringType), Seq.empty)

    womValueToMetadataEvents(MetadataKey(workflowId, None, "root"), womArray).toList should contain theSameElementsInOrderAs List(
      MetadataEvent(MetadataKey(workflowId, None, "root[0]"), MetadataValue("Hello")),
      MetadataEvent(MetadataKey(workflowId, None, "root[1]"), MetadataValue("world!"))
    )

    womValueToMetadataEvents(MetadataKey(workflowId, None, "root"), emptyWomArray).toList should contain theSameElementsAs List(
      MetadataEvent.empty(MetadataKey(workflowId, None, "root[]"))
    )
  }

  it should "convert a WomMap to MetadataEvents" in {
    import MetadataService._
    val workflowId = WorkflowId.randomId()
    val womMap = WomMap(WomMapType(WomStringType, WomStringType), Map(
      WomString("Hello") -> WomString("world!"),
      WomString("Goodbye") -> WomString("world!")
    ))
    val emptyWomMap = WomMap(WomMapType(WomStringType, WomStringType), Map.empty)

    womValueToMetadataEvents(MetadataKey(workflowId, None, "root"), womMap).toList should contain theSameElementsInOrderAs List(
      MetadataEvent(MetadataKey(workflowId, None, "root:Hello"), MetadataValue("world!")),
      MetadataEvent(MetadataKey(workflowId, None, "root:Goodbye"), MetadataValue("world!"))
    )

    womValueToMetadataEvents(MetadataKey(workflowId, None, "root"), emptyWomMap).toList should contain theSameElementsAs List(
      MetadataEvent.empty(MetadataKey(workflowId, None, "root"))
    )
  }

  it should "convert a WomMaybePopulatedFile to MetadataEvents" in {
    import MetadataService._
    val workflowId = WorkflowId.randomId()
    // Leave some fields out on purpose to verify that they still make it (empty) to the metadata
    val womMaybePopulatedFileInner = WomMaybePopulatedFile(
      valueOption = Option("innerFileValue"),
      contentsOption = Option("innerContents")
    )

    val womMaybePopulatedFileOuter = WomMaybePopulatedFile(
      valueOption = Option("fileValue"),
      checksumOption = Option("checksum"),
      sizeOption = Option(5L),
      formatOption = Option("format"),
      contentsOption = Option("contents"),
      secondaryFiles = List(
        womMaybePopulatedFileInner,
        WomSingleFile("outerSingleFile")
      )
    )

    womValueToMetadataEvents(MetadataKey(workflowId, None, "root"), womMaybePopulatedFileOuter).toList should contain theSameElementsInOrderAs List(
      MetadataEvent(MetadataKey(workflowId, None, "root:value"), MetadataValue("fileValue")),
      MetadataEvent(MetadataKey(workflowId, None, "root:checksum"), MetadataValue("checksum")),
      MetadataEvent(MetadataKey(workflowId, None, "root:size"), MetadataValue(5L)),
      MetadataEvent(MetadataKey(workflowId, None, "root:format"), MetadataValue("format")),
      MetadataEvent(MetadataKey(workflowId, None, "root:contents"), MetadataValue("contents")),
      MetadataEvent(MetadataKey(workflowId, None, "root:secondaryFiles[0]:value"), MetadataValue("innerFileValue")),
      MetadataEvent(MetadataKey(workflowId, None, "root:secondaryFiles[0]:checksum"), None),
      MetadataEvent(MetadataKey(workflowId, None, "root:secondaryFiles[0]:size"), None),
      MetadataEvent(MetadataKey(workflowId, None, "root:secondaryFiles[0]:format"), None),
      MetadataEvent(MetadataKey(workflowId, None, "root:secondaryFiles[0]:contents"), MetadataValue("innerContents")),
      MetadataEvent(MetadataKey(workflowId, None, "root:secondaryFiles[0]:secondaryFiles[]"), None),
      MetadataEvent(MetadataKey(workflowId, None, "root:secondaryFiles[1]"), MetadataValue("outerSingleFile"))
    )
  }

  it should "convert a WomMaybeListedDirectory to MetadataEvents" in {
    import MetadataService._
    val workflowId = WorkflowId.randomId()
    // Leave some fields out on purpose to verify that they still make it (empty) to the metadata
    val womMaybeListedDirectoryInner = WomMaybeListedDirectory(
      valueOption = Option("innerDirectoryValue"),
      listingOption = None
    )

    val womMaybeListedDirectoryOuter = WomMaybeListedDirectory(
      valueOption = Option("directoryValue"),
      listingOption = Option(
        List(
          womMaybeListedDirectoryInner,
          WomSingleFile("outerSingleFile")
        )
      )
    )

    womValueToMetadataEvents(MetadataKey(workflowId, None, "root"), womMaybeListedDirectoryOuter).toList should contain theSameElementsInOrderAs List(
      MetadataEvent(MetadataKey(workflowId, None, "root:value"), MetadataValue("directoryValue")),
      MetadataEvent(MetadataKey(workflowId, None, "root:listing[0]:value"), MetadataValue("innerDirectoryValue")),
      MetadataEvent(MetadataKey(workflowId, None, "root:listing[0]:listing[]"), None),
      MetadataEvent(MetadataKey(workflowId, None, "root:listing[1]"), MetadataValue("outerSingleFile"))
    )
  }

  it should "convert a primitive WomValue to MetadataEvents" in {
    import MetadataService._
    val workflowId = WorkflowId.randomId()

    val values = Table(
      ("womValue", "metadataValue"),
      (WomString("hi"), MetadataValue("hi", MetadataString)),
      (WomInteger(1), MetadataValue("1", MetadataInt)),
      (WomFloat(1F), MetadataValue("1.0", MetadataNumber)),
      (WomBoolean(true), MetadataValue("true", MetadataBoolean))
    )

    forAll(values) { (womValue, metadataValue) =>
      womValueToMetadataEvents(MetadataKey(workflowId, None, "root"), womValue).toList should contain theSameElementsAs List(
        MetadataEvent(MetadataKey(workflowId, None, "root"), metadataValue)
      )
    }
  }

  // For the metadata event tests!
  val failureMessageRegex = "([^\\[]*)\\[([0-9]+)\\](.*)\\:message".r
  val pathToFailures = "path:to:failures"

  it should "convert an exception into a failure event with an empty causedBy block" in {
    import MetadataService.throwableToMetadataEvents

    val workflowId = WorkflowId(UUID.randomUUID())
    val mdkey = MetadataKey(workflowId, None, pathToFailures)

    val tMsg = "The Oscars suck!"
    val t = new RuntimeException(tMsg)

    val events = throwableToMetadataEvents(mdkey, t)
    events.size should be(2)
    val (keyPrefix, causedBys, failureIndex) = validateExceptionMessage(events.head, workflowId, tMsg)
    keyPrefix should be(pathToFailures)
    causedBys should be("")
    events(1).key.key should be(s"$keyPrefix[$failureIndex]:causedBy[]")
    events(1).key.workflowId should be(workflowId)
    events(1).value should be(None)

  }

  it should "convert nested exceptions into a sequence of failure events" in {
    import MetadataService.throwableToMetadataEvents

    val workflowId = WorkflowId(UUID.randomUUID())
    val mdkey = MetadataKey(workflowId, None, pathToFailures)

    val innerCauseMsg = "Envelope malfunctions"
    val innerCause = new RuntimeException(innerCauseMsg)

    val causeMsg = "Wrong recipients"
    val cause = new RuntimeException(causeMsg, innerCause)

    val tMsg = "The Oscars suck!"
    val t = new RuntimeException(tMsg, cause)

    val events = throwableToMetadataEvents(mdkey, t)
    events.size should be(4)

    val (outerPrefix, outerCausedBys, outerFailureId) = validateExceptionMessage(events.head, workflowId, tMsg)
    val (cause1Prefix, cause1CausedBys, cause1FailureId) = validateExceptionMessage(events(1), workflowId, causeMsg)
    val (cause2Prefix, cause2CausedBys, cause2FailureId) = validateExceptionMessage(events(2), workflowId, innerCauseMsg)
    events(3).key.key should be(s"$cause2Prefix[$cause2FailureId]$cause2CausedBys:causedBy[]")

    outerPrefix should be(pathToFailures)
    cause1Prefix should be(pathToFailures)
    cause2Prefix should be(pathToFailures)
    outerCausedBys should be("")
    cause1CausedBys should be(":causedBy[0]")
    cause2CausedBys should be(":causedBy[0]:causedBy[0]")
    cause1FailureId should be(outerFailureId)
    cause2FailureId should be(cause1FailureId)
  }

  it should "convert aggregated exceptions into a sequence of failure events" in {
    import MetadataService.throwableToMetadataEvents

    val workflowId = WorkflowId(UUID.randomUUID())
    val mdkey = MetadataKey(workflowId, None, "path:to:failures")

    val innerCauseMsg = "Envelope malfunctions"
    val innerCause = new RuntimeException(innerCauseMsg)

    val cause1Msg = "Wrong recipients"
    val cause1 = new RuntimeException(cause1Msg)
    val cause2Msg = "Self congratulation"
    val cause2 = new RuntimeException(cause2Msg, innerCause)
    val cause3Msg = "The Globes are better anyway"
    val cause3 = new RuntimeException(cause3Msg)

    val causeContext = "Compound Entertainment Failure"
    val cause = new AggregatedException(causeContext, List(cause1, cause2, cause3))

    val tMsg = "The Oscars suck!"
    val t = new RuntimeException(tMsg, cause)

    val events = throwableToMetadataEvents(mdkey, t)
    events.size should be(9)

    // Outer runtime exception:
    val (runtimePrefix, runtimeCausedBys, runtimeFailureId) = validateExceptionMessage(events.head, workflowId, tMsg)
    runtimePrefix should be(pathToFailures)
    runtimeCausedBys should be("")

    // Aggregate exception:
    val (aggregatePrefix, aggregateCausedBys, aggregateFailureId) = validateExceptionMessage(events(1), workflowId, causeContext)
    aggregatePrefix should be(pathToFailures)
    aggregateCausedBys should be(":causedBy[0]")
    aggregateFailureId should be(runtimeFailureId)

    // cause1, caused by []
    val (cause1Prefix, cause1CausedBys, cause1FailureId) = validateExceptionMessage(events(2), workflowId, cause1Msg)
    cause1Prefix should be(pathToFailures)
    cause1CausedBys should be(":causedBy[0]:causedBy[0]")
    cause1FailureId should be(runtimeFailureId)
    events(3).key.key should be(s"$cause1Prefix[$runtimeFailureId]$cause1CausedBys:causedBy[]")

    // cause2, caused by innerCause caused by []
    val (cause2Prefix, cause2CausedBys, cause2FailureId) = validateExceptionMessage(events(4), workflowId, cause2Msg)
    val (innerCausePrefix, innerCauseCausedBys, innerCauseFailureIds) = validateExceptionMessage(events(5), workflowId, innerCauseMsg)
    cause2Prefix should be(pathToFailures)
    cause2CausedBys should be(":causedBy[0]:causedBy[1]")
    cause2FailureId should be(runtimeFailureId)
    innerCausePrefix should be(pathToFailures)
    innerCauseCausedBys should be(":causedBy[0]:causedBy[1]:causedBy[0]")
    innerCauseFailureIds should be(runtimeFailureId)
    events(6).key.key should be(s"$innerCausePrefix[$runtimeFailureId]$innerCauseCausedBys:causedBy[]")

    // cause3, caused by []
    val (cause3Prefix, cause3CausedBys, cause3FailureId) = validateExceptionMessage(events(7), workflowId, cause3Msg)
    cause3Prefix should be(pathToFailures)
    cause3CausedBys should be(":causedBy[0]:causedBy[2]")
    cause3FailureId should be(runtimeFailureId)
    events(8).key.key should be(s"$cause3Prefix[$cause3FailureId]$cause3CausedBys:causedBy[]")
  }

  def validateExceptionMessage(event: MetadataEvent, workflowId: WorkflowId, message: String) = event match {
    case MetadataEvent(k, Some(MetadataValue(v, _)), _) =>
      k.workflowId should be(workflowId)
      v should be(message)

      // Return the ID so that we can check for uniqueness later:
      k.key match {
        case failureMessageRegex(prefix, failureIndex, causedBys) => (prefix, causedBys, failureIndex)
        case _ => fail("Unexpected failure key format: " + k.key)
      }
    case _ => fail("throwableToMetadataEvents generated a metadata event without a metadata value! Bad throwableToMetadataEvents! Very bad!")
  }
}
