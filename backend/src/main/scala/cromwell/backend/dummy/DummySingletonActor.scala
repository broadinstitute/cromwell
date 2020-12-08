package cromwell.backend.dummy

import java.io.File
import java.time.OffsetDateTime
import java.util.UUID

import akka.actor.Actor
import com.typesafe.scalalogging.StrictLogging
import cromwell.backend.dummy.DummySingletonActor._

import scala.concurrent.ExecutionContext
import scala.concurrent.duration._

final class DummySingletonActor() extends Actor with StrictLogging {

  implicit val ec: ExecutionContext = context.dispatcher
  var count: Int = 0

  var countHistory: Vector[(OffsetDateTime, Int)] = Vector.empty

  override def receive: Receive = {
    case PlusOne => count = count + 1
    case MinusOne => count = count - 1
    case PrintCount =>
      if(countHistory.lastOption.exists(_._2 != count)) {
        countHistory = countHistory :+ (OffsetDateTime.now() -> count)
        logger.info("The current count is now: " + count)
        if (count == 0) {
          outputCountHistory()
          countHistory = Vector.empty
        }
      } else {
        countHistory = countHistory :+ (OffsetDateTime.now() -> count)
      }

  }

  private def outputCountHistory() = {
    import java.io.BufferedWriter
    import java.io.FileOutputStream
    import java.io.OutputStreamWriter
    val fout = new File(s"timestamps-${UUID.randomUUID().toString}.tsv")
    val fos = new FileOutputStream(fout)

    val bw = new BufferedWriter(new OutputStreamWriter(fos))

    for ((timestamp, count) <- countHistory) {
      bw.write(s"$timestamp\t$count")
      bw.newLine()
    }
    bw.flush()
    bw.close()
  }

  context.system.scheduler.schedule(10.seconds, 1.second) { self ! PrintCount }
}

object DummySingletonActor {
  case object PlusOne
  case object MinusOne
  case object PrintCount
}

