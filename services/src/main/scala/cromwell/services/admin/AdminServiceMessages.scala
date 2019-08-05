package cromwell.services.admin

import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

object AdminServiceMessages {

  final val AdminServiceName = "Admin"

  sealed trait AdminServiceMessage extends ServiceRegistryMessage {
    override def serviceName: String = AdminServiceName
  }

  case object ListSubmissionsRequest extends AdminServiceMessage
  case object PauseSubmissionRequest extends AdminServiceMessage


  sealed trait AdminResult extends AdminServiceMessage
  case object PauseSubmissionSuccess extends AdminResult
}

