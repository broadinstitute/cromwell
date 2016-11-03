package cromwell.backend.impl.jes

import java.time.OffsetDateTime
import java.util

import com.google.api.client.util.ArrayMap
import com.google.api.services.genomics.model.Operation
import cromwell.backend.impl.jes.Run.GlobInfo
import cromwell.core.ExecutionEvent
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.{Mockito => MockitoTrait}

import scala.collection.JavaConverters._

class RunSpec extends FlatSpec with Matchers with MockitoTrait {
  "JES Run" should "parse events from Operation metadata" in {
    val op: Operation = new Operation()

    // An execution event
    val executionEvent: ArrayMap[String, String] = ArrayMap.create(2)
    executionEvent.add("description", "start")
    executionEvent.add("startTime", "2015-12-05T00:00:01+00:00")

    // An event that we should ignore
    val eventToIgnore: ArrayMap[String, String] = ArrayMap.create(2)
    eventToIgnore.add("description", "blah")
    eventToIgnore.add("startTime", "2015-12-05T00:01:00+00:00")

    // An event we should pick out as a glob
    val globEvent: ArrayMap[String, String] = ArrayMap.create(2)
    globEvent.add("description", "copied 100 file(s) to \"gs://some_path/glob-7bb0c33ac658a900e2c3804726fc1d2a/\"")
    globEvent.add("startTime", "2015-12-05T00:01:00+00:00") // Completely ignored, but still needed

    val events = new util.ArrayList(Seq(executionEvent, eventToIgnore, globEvent).asJava)

    val metadata: Map[String, AnyRef] = Map(
      "createTime" -> "2015-12-05T00:00:00+00:00",
      "startTime" -> "2015-12-05T00:00:01+00:00",
      "endTime" -> "2015-12-05T11:00:00+00:00",
      "events" -> events
    )

    op.setMetadata(metadata.asJava)

    val list = Run.getEventList(op)
    list.events should contain theSameElementsAs List(
      ExecutionEvent("waiting for quota", OffsetDateTime.parse("2015-12-05T00:00:00+00:00")),
      ExecutionEvent("initializing VM", OffsetDateTime.parse("2015-12-05T00:00:01+00:00")),
      ExecutionEvent("start", OffsetDateTime.parse("2015-12-05T00:00:01+00:00")),
      ExecutionEvent("cromwell poll interval", OffsetDateTime.parse("2015-12-05T11:00:00+00:00"))
    )

    list.globs.size should be(1)
    list.globs.head should be(GlobInfo("gs://some_path/glob-7bb0c33ac658a900e2c3804726fc1d2a/", 100))
  }
}
