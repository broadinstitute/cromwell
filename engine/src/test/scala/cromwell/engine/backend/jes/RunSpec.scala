package cromwell.engine.backend.jes

import java.time.OffsetDateTime
import java.util

import com.google.api.client.googleapis.testing.auth.oauth2.MockGoogleCredential
import com.google.api.client.util.ArrayMap
import com.google.api.services.genomics.Genomics
import com.google.api.services.genomics.model.Operation
import cromwell.backend.impl.jes.{EventStartTime, Run}
import org.scalatest.{FlatSpec, Matchers}
import org.specs2.mock.{Mockito => MockitoTrait}

import scala.collection.JavaConverters._

class RunSpec extends FlatSpec with Matchers with MockitoTrait {
  "JES Run" should "parse events from Operation metadata" in {
    val op: Operation = new Operation()

    val event1: ArrayMap[String, String] = ArrayMap.create(2)
    event1.add("description", "blah")
    event1.add("startTime", "2015-12-05T00:00:01+00:00")

    val event2: ArrayMap[String, String] = ArrayMap.create(2)
    event2.add("description", "blah2")
    event2.add("startTime", "2015-12-05T00:01:00+00:00")

    val events = new util.ArrayList(Seq(event1, event2).asJava)

    val metadata: Map[String, AnyRef] = Map(
      "createTime" -> "2015-12-05T00:00:00+00:00",
      "startTime" -> "2015-12-05T00:00:01+00:00",
      "endTime" -> "2015-12-05T11:00:00+00:00",
      "events" -> events
    )

    op.setMetadata(metadata.asJava)

    val mockedCredentials = new MockGoogleCredential.Builder().build()
    val genomics = new Genomics(mockedCredentials.getTransport, mockedCredentials.getJsonFactory, mockedCredentials)
    val run = new Run("runId", genomics)
    val list = run.getEventList(op)
    list should contain theSameElementsAs List(
      EventStartTime("waiting for quota", OffsetDateTime.parse("2015-12-05T00:00:00+00:00")),
      EventStartTime("initializing VM", OffsetDateTime.parse("2015-12-05T00:00:01+00:00")),
      EventStartTime("blah", OffsetDateTime.parse("2015-12-05T00:00:01+00:00")),
      EventStartTime("blah2", OffsetDateTime.parse("2015-12-05T00:01:00+00:00")),
      EventStartTime("cromwell poll interval", OffsetDateTime.parse("2015-12-05T11:00:00+00:00"))
    )

  }
}
