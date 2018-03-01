package cromwell.services.loadcontroller

import cromwell.services.ServiceRegistryActor.{ListenToMessage, ServiceRegistryMessage}
import cromwell.services.loadcontroller.impl.LoadControllerServiceActor.LoadControllerServiceName

object LoadControllerService {
  sealed trait LoadControllerMessage extends ServiceRegistryMessage {
    def serviceName = LoadControllerServiceName
  }
  sealed trait LoadLevel { def level: Int }
  case object NormalLoad extends LoadLevel { val level = 0 }
  // Start with a single abnormal load level for now and we can add more if we want to be more granular
  case object HighLoad extends LoadLevel { val level = 1 }

  implicit val loadLevelOrdering: Ordering[LoadLevel] = Ordering.by[LoadLevel, Int](_.level)
  case object ListenToLoadController extends LoadControllerMessage with ListenToMessage
}
