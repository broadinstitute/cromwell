package cromwell.server

import akka.actor.ActorSystem
import akka.event.EventStream
import akka.event.slf4j.Slf4jLoggingFilter

class CromwellAkkaLogFilter(settings: ActorSystem.Settings, eventStream: EventStream)
    extends Slf4jLoggingFilter(settings, eventStream) {
  override def isErrorEnabled(logClass: Class[_], logSource: String) =
    /*
     * This might filter out too much but it's the finest granularity we have here
     * The goal is to not log the
     * "Outgoing request stream error akka.stream.AbruptTerminationException:
     *  Processor actor [Actor[akka://cromwell-system/user/StreamSupervisor-1/flow-6-0-mergePreferred#1200284127]] terminated abruptly"
     * type of message
     *
     * See https://github.com/akka/akka-http/issues/907 and https://github.com/akka/akka/issues/18747
     */
    super.isErrorEnabled(logClass, logSource) && !(logSource.startsWith(
      "akka.actor.ActorSystemImpl"
    ) && CromwellShutdown.shutdownInProgress())
}
