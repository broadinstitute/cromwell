package cromwell.engine.backend.jes

import java.util

import com.google.api.client.util.ArrayMap
import com.google.api.services.genomics.model.Operation
import org.joda.time.DateTime
import org.scalatest.{Matchers, FlatSpec}
import scala.collection.JavaConverters._


class RunSpec extends FlatSpec with Matchers {
  "JES Run" should "parse events from Operation metadata" in {
    val op: Operation = new Operation()

    val event1: ArrayMap[String, String] = ArrayMap.create(2)
    event1.add("description", "pulling-image")
    event1.add("startTime", "2015-12-05T00:00:01+00:00")

    val event2: ArrayMap[String, String] = ArrayMap.create(2)
    event2.add("description", "start")
    event2.add("startTime", "2015-12-05T00:01:00+00:00")

    val events = new util.ArrayList(Seq(event1, event2).asJava)

    val metadata: Map[String, AnyRef] = Map(
      "createTime" -> "2015-12-05T00:00:00+00:00",
      "startTime" -> "2015-12-05T00:00:01+00:00",
      "endTime" -> "2015-12-05T11:00:00+00:00",
      "events" -> events
    )

    op.setMetadata(metadata.asJava)

    Run.getEventList(op) should have size 5
    Run.getEventList(op) filter { x => x.description == "pulling-image" } foreach { x =>
      x.startTime.getMillis should be (new DateTime("2015-12-05T00:00:01.000Z").getMillis)
      x.endTime.getMillis should be (new DateTime("2015-12-05T00:01:00.000Z").getMillis)
    }
  }

  "JES Run" should "not parse unwated events from Operation metadata" in {
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

    Run.getEventList(op) should have size 3
    Run.getEventList(op) filter { x => x.description == "blah" } foreach { x =>
      x.startTime.getMillis should be (None)
      x.endTime.getMillis should be (None)
    }
  }
}
