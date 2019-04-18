package cromwell.backend.google.pipelines.v1alpha2.api

import java.time.OffsetDateTime
import java.util

import com.google.api.client.util.ArrayMap
import com.google.api.services.genomics.model.Operation
import cromwell.backend.google.pipelines.common.api.RunStatus.Success
import cromwell.backend.google.pipelines.v1alpha2.api.request.GetRequestHandler
import cromwell.core.ExecutionEvent
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.{Mockito => MockitoTrait}

import scala.collection.JavaConverters._

class StatusResponseParsingSpec extends FlatSpec with Matchers with MockitoTrait {
  behavior of "Pipelines Api status parsing"

  it should "parse events from Operation metadata" in {
    val op: Operation = new Operation()
    op.setMetadata(eventsMetadata.asJava)

    val list = GetRequestHandler.getEventList(op)
    list should contain theSameElementsAs eventsExpected
  }

  it should "require operation be non-null" in {
    val exception = intercept[RuntimeException](GetRequestHandler.interpretOperationStatus(null))
    exception.getMessage should be("requirement failed: Operation must not be null.")
  }

  it should "catch and wrap null pointer exceptions in an empty operation" in {
    val op = new Operation()

    val exception = intercept[RuntimeException](GetRequestHandler.interpretOperationStatus(op))
    exception.getMessage should startWith("Caught NPE while processing operation null")
  }

  it should "catch and wrap null pointer exceptions in a name only operation" in {
    val op = new Operation()
    op.setName("my/customName")

    val exception = intercept[RuntimeException](GetRequestHandler.interpretOperationStatus(op))
    exception.getMessage should startWith ("Caught NPE while processing operation my/customName")
  }

  it should "parse an operation without machine information" in {
    val op = new Operation()
    op.setName("my/customName")
    op.setDone(true)
    op.setMetadata(eventsMetadata.asJava)

    val runStatus = GetRequestHandler.interpretOperationStatus(op)

    runStatus should be(a[Success])

    val success = runStatus.asInstanceOf[Success]
    success.instanceName should be(None)
    success.machineType should be(None)
    success.zone should be(None)
    success.eventList should contain theSameElementsAs eventsExpected
  }

  it should "parse an operation with empty runtime metadata" in {
    val op = new Operation()
    op.setName("my/customName")
    op.setDone(true)

    val runtimeMetadata = ArrayMap.create[String, Object]()

    val metadata = eventsMetadata ++ Map("runtimeMetadata" -> runtimeMetadata)
    op.setMetadata(metadata.asJava)

    val runStatus = GetRequestHandler.interpretOperationStatus(op)
    runStatus should be(a[Success])

    val success = runStatus.asInstanceOf[Success]
    success.instanceName should be(None)
    success.machineType should be(None)
    success.zone should be(None)
    success.eventList should contain theSameElementsAs eventsExpected
  }

  it should "parse an operation with empty compute engine information" in {
    val op = new Operation()
    op.setName("my/customName")
    op.setDone(true)

    val computeEngine = ArrayMap.create[String, String]()

    val runtimeMetadata = ArrayMap.create[String, Object]()
    runtimeMetadata.add("computeEngine", computeEngine)

    val metadata = eventsMetadata ++ Map("runtimeMetadata" -> runtimeMetadata)
    op.setMetadata(metadata.asJava)

    val runStatus = GetRequestHandler.interpretOperationStatus(op)
    runStatus should be(a[Success])

    val success = runStatus.asInstanceOf[Success]
    success.instanceName should be(None)
    success.machineType should be(None)
    success.zone should be(None)
    success.eventList should contain theSameElementsAs eventsExpected
  }

  it should "parse an operation with partially filled compute engine" in {
    val op = new Operation()
    op.setName("my/customName")
    op.setDone(true)

    val computeEngine = ArrayMap.create[String, String]()
    computeEngine.add("zone", "us-central1-b")

    val runtimeMetadata = ArrayMap.create[String, Object]()
    runtimeMetadata.add("computeEngine", computeEngine)

    val metadata = eventsMetadata ++ Map("runtimeMetadata" -> runtimeMetadata)
    op.setMetadata(metadata.asJava)

    val runStatus = GetRequestHandler.interpretOperationStatus(op)

    runStatus should be(a[Success])

    val success = runStatus.asInstanceOf[Success]
    success.instanceName should be(None)
    success.machineType should be(None)
    success.zone should be(Option("us-central1-b"))
    success.eventList should contain theSameElementsAs eventsExpected
  }

  it should "parse an operation with filled compute engine" in {
    val op = new Operation()
    op.setName("my/customName")
    op.setDone(true)

    val computeEngine = ArrayMap.create[String, String]()
    computeEngine.add("instanceName", "ggp-12345678901234567890")
    computeEngine.add("machineType", "us-central1-b/g1-small")
    computeEngine.add("zone", "us-central1-b")

    val runtimeMetadata = ArrayMap.create[String, Object]()
    runtimeMetadata.add("computeEngine", computeEngine)

    val metadata = eventsMetadata ++ Map("runtimeMetadata" -> runtimeMetadata)
    op.setMetadata(metadata.asJava)

    val runStatus = GetRequestHandler.interpretOperationStatus(op)

    runStatus should be(a[Success])

    val success = runStatus.asInstanceOf[Success]
    success.instanceName should be(Option("ggp-12345678901234567890"))
    success.machineType should be(Option("us-central1-b/g1-small"))
    success.zone should be(Option("us-central1-b"))
    success.eventList should contain theSameElementsAs eventsExpected
  }

  private lazy val eventsMetadata: Map[String, AnyRef] = {
    val event1: ArrayMap[String, String] = ArrayMap.create(2)
    event1.add("description", "start")
    event1.add("startTime", "2015-12-05T00:00:01+00:00")

    val event2: ArrayMap[String, String] = ArrayMap.create(2)
    event2.add("description", "blah")
    event2.add("startTime", "2015-12-05T00:01:00+00:00")

    val events = new util.ArrayList(Seq(event1, event2).asJava)

    Map(
      "createTime" -> "2015-12-05T00:00:00+00:00",
      "startTime" -> "2015-12-05T00:00:01+00:00",
      "endTime" -> "2015-12-05T11:00:00+00:00",
      "events" -> events
    )
  }

  private lazy val eventsExpected =
    List(
      ExecutionEvent("waiting for quota", OffsetDateTime.parse("2015-12-05T00:00:00+00:00")),
      ExecutionEvent("initializing VM", OffsetDateTime.parse("2015-12-05T00:00:01+00:00")),
      ExecutionEvent("start", OffsetDateTime.parse("2015-12-05T00:00:01+00:00")),
      ExecutionEvent("cromwell poll interval", OffsetDateTime.parse("2015-12-05T11:00:00+00:00"))
    )
}
