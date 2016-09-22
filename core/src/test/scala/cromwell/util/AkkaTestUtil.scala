package cromwell.util

import akka.actor.{Actor, ActorLogging, Props}
import akka.testkit.TestProbe

object AkkaTestUtil {

  implicit class EnhancedTestProbe(probe: TestProbe) {
    def props = Props(new Actor with ActorLogging {
      def receive = {
        case outbound if sender == probe.ref =>
          val msg = "Unexpected outbound message from Probe. You're doing something wrong!"
          log.error(msg)
          throw new RuntimeException(msg)
        case inbound => probe.ref forward inbound
      }
    })
  }
}
