package cromwell.services.cost

import akka.actor.{Actor, ActorRef}
import com.google.cloud.billing.v1._
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.util.StringUtil.EnhancedToStringable
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

case class CostCatalogMessage() extends ServiceRegistryMessage {
  override val serviceName = "CostCatalogService"
}

class CostCatalogService(serviceConfig: Config, globalConfig: Config, serviceRegistry: ActorRef) extends Actor with LazyLogging{
  override def receive: Receive = {
    case CostCatalogMessage() => {
      fetchPublicCostCatalog()
      logger.error(
        s"Got message"
      )
    }
    case ShutdownCommand => context stop self
    case other =>
      logger.error(
        s"Programmer Error: Unexpected message ${other.toPrettyElidedString(1000)} received by ${this.self.path.name}."
      )
  }
  def serviceRegistryActor: ActorRef = serviceRegistry

  def fetchPublicCostCatalog(): Unit = {
    val cloudCatalogClient = CloudCatalogClient.create
    val request = ListSkusRequest.newBuilder.setParent(ServiceName.of("[SERVICE]").toString).setCurrencyCode("currencyCode1004773790").setPageSize(883849137).setPageToken("pageToken873572522").build
    val response = cloudCatalogClient.listSkus(request)
    response.iterateAll().forEach(println)
  }

}
