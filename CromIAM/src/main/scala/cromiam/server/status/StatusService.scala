package cromiam.server.status

import akka.actor.ActorSystem
import akka.pattern.ask
import akka.util.Timeout
import org.broadinstitute.dsde.workbench.util.health.HealthMonitor.GetCurrentStatus
import org.broadinstitute.dsde.workbench.util.health.{HealthMonitor, StatusCheckResponse, SubsystemStatus}
import org.broadinstitute.dsde.workbench.util.health.Subsystems._

import scala.concurrent.{ExecutionContext, Future}
import scala.concurrent.duration._

/**
  * Sets up and provides access to the health monitor for CromIAM
  */
class StatusService(checkStatus: () => Map[Subsystem, Future[SubsystemStatus]],
                    initialDelay: FiniteDuration = Duration.Zero,
                    pollInterval: FiniteDuration = 1.minute)(implicit system: ActorSystem, executionContext: ExecutionContext) {
  implicit val askTimeout = Timeout(5.seconds)

  private val healthMonitor = system.actorOf(HealthMonitor.props(Set(Cromwell, Sam))(checkStatus), "HealthMonitorActor")

  system.scheduler.scheduleWithFixedDelay(initialDelay, pollInterval, healthMonitor, HealthMonitor.CheckAll)

  def status(): Future[StatusCheckResponse] = healthMonitor.ask(GetCurrentStatus).asInstanceOf[Future[StatusCheckResponse]]
}

