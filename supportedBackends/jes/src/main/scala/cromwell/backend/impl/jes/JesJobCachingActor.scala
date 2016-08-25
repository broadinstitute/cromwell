package cromwell.backend.impl.jes

import akka.actor.{Props, ActorRef}
import cromwell.backend.BackendJobCachingActor.BackendJobExecutionResponse
import cromwell.backend.{BackendJobCachingActor, BackendJobDescriptor}
import org.slf4j.LoggerFactory
import scala.concurrent.Future

object JesJobCachingActor {
  val logger = LoggerFactory.getLogger("JesBackend")

  def props(jobDescriptor: BackendJobDescriptor,
            jesConfiguration: JesConfiguration,
            initializationData: JesBackendInitializationData,
            serviceRegistryActor: ActorRef): Props = {
    Props(new JesJobCachingActor(jobDescriptor, jesConfiguration, initializationData, serviceRegistryActor))
  }
}

case class JesJobCachingActor(override val jobDescriptor: BackendJobDescriptor,
                                jesConfiguration: JesConfiguration,
                                initializationData: JesBackendInitializationData,
                                serviceRegistryActor: ActorRef)
  extends BackendJobCachingActor {

  override def abort = ???

  override def copyCachedOutputs(): Future[BackendJobExecutionResponse]
}