package cromwell.services.admin

import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage

object AdminServiceMessages {

  final val AdminServiceName = "Admin"

  sealed trait AdminServiceMessage extends ServiceRegistryMessage {
    override def serviceName: String = AdminServiceName
  }

  case object ListSubmissionsRequest extends AdminServiceMessage
  case object PauseSubmissionRequest extends AdminServiceMessage

  sealed trait PauseSubmissionResult extends AdminServiceMessage
  case object PauseSubmissionSuccess extends PauseSubmissionResult
  case class PauseSubmissionFailure(reason: String) extends PauseSubmissionResult

//  sealed trait ListSubmissionsResult extends AdminServiceMessage
//  case class ListSubmissionsSuccess(submissions: List[ListSubmissionData]) extends ListSubmissionsResult
//  case class ListSubmissionsFailure(reason: Throwable) extends ListSubmissionsResult
}
