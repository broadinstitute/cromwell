package cromwell.services.cost

import akka.actor.{Actor, ActorRef}
import com.google.cloud.billing.v1._
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.util.StringUtil.EnhancedToStringable
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import java.io.{FileOutputStream, ObjectOutputStream}
import scala.jdk.CollectionConverters.IterableHasAsScala
case class CostCatalogMessage() extends ServiceRegistryMessage {
  override val serviceName = "CostCatalogService"
}

class CostCatalogService(serviceConfig: Config, globalConfig: Config, serviceRegistry: ActorRef)
    extends Actor
    with LazyLogging {
  private val COMPUTE_ENGINE_SERVICE_NAME =
    "services/6F81-5844-456A" // Can be gleaned by using googleClient.listServices
  private val DEFAULT_CURRENCY_CODE =
    "USD" // ISO 4217 https://developers.google.com/adsense/management/appendix/currencies
  private lazy val googleClient = CloudCatalogClient.create
  private lazy val costCatalog: Map[CostCatalogKey, CostCatalogValue] = processCostCatalog(makeInitialWebRequest().iterateAll().asScala)

  def getCatalog: Map[CostCatalogKey, CostCatalogValue] = costCatalog
  def saveResponseToFileForTesting(filepath: String) = {
    val dataToSave: List[Sku] = makeInitialWebRequest().iterateAll().asScala.toList
    val fileOutputStream = new FileOutputStream(filepath)
    val objectOutputStream = new ObjectOutputStream(fileOutputStream)
    objectOutputStream.writeObject(dataToSave)
    objectOutputStream.flush()
  }

  override def receive: Receive = {
    case ShutdownCommand => context stop self
    case other =>
      logger.error(
        s"Programmer Error: Unexpected message ${other.toPrettyElidedString(1000)} received by ${this.self.path.name}."
      )
  }
  def serviceRegistryActor: ActorRef = serviceRegistry

  /**
   * Makes a web request to google that lists all SKUs. The request is paginated, so the response will only include
   * the first page. In order to get to all SKUs via iterateAll(), additional web requests will be made.
   */
  private def makeInitialWebRequest() : CloudCatalogClient.ListSkusPagedResponse = {
    val request = ListSkusRequest
      .newBuilder()
      .setParent(COMPUTE_ENGINE_SERVICE_NAME)
      .setCurrencyCode(DEFAULT_CURRENCY_CODE)
      .build()
    googleClient.listSkus(request)
  }

  private def processCostCatalog(skus: Iterable[Sku]): Map[CostCatalogKey, CostCatalogValue] = {
    // TODO: Account for key collisions (same key can be in multiple regions)
    // TODO: reduce memory footprint of returned map  (don't store entire SKU object)
    skus.foldLeft(Map.empty[CostCatalogKey, CostCatalogValue]) { case (acc, sku) =>
      acc + CostCatalogUtils.convertSkuToKeyValuePair(sku)
    }
  }
}
