package cromiam.server.status

import akka.actor.ActorSystem
import org.broadinstitute.dsde.workbench.util.health.{StatusCheckResponse, SubsystemStatus}
import org.broadinstitute.dsde.workbench.util.health.Subsystems.{Cromwell, Sam, Subsystem}

import scala.concurrent.{ExecutionContext, Future}

class MockStatusService(checkStatus: () => Map[Subsystem, Future[SubsystemStatus]])(implicit
  system: ActorSystem,
  executionContext: ExecutionContext
) extends StatusService(checkStatus)(system, executionContext) {

  override def status(): Future[StatusCheckResponse] = {
    val subsystemStatus: SubsystemStatus = SubsystemStatus(ok = true, None)
    val subsystems: Map[Subsystem, SubsystemStatus] = Map(Cromwell -> subsystemStatus, Sam -> subsystemStatus)

    Future.successful(StatusCheckResponse(ok = true, subsystems))
  }
}
