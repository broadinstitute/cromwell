package cromwell.services.cost

import akka.actor.{Actor, ActorRef}
import com.google.cloud.billing.v1._
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.util.StringUtil.EnhancedToStringable
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import java.io.{FileOutputStream, ObjectOutputStream}
import java.time.{Duration, LocalDateTime}
import scala.jdk.CollectionConverters.IterableHasAsScala
import java.time.temporal.ChronoUnit.HOURS
case class CostCatalogMessage() extends ServiceRegistryMessage {
  override val serviceName = "CostCatalogService"
}

case class ExpiringCostCatalog(catalog: Map[CostCatalogKey, CostCatalogValue], fetchTime: LocalDateTime)

/**
 * This actor handles the fetching and caching of the public Google Cost Catalog.
 * This enables us to look up the SKUs necessary for cost calculations.
 */
class CostCatalogService(serviceConfig: Config,
                         globalConfig: Config,
                         serviceRegistry: ActorRef,
                         mockTestData: Option[List[Sku]] = Option.empty,
                         maxCatalogLifetime: Duration = Duration.of(24, HOURS)
) extends Actor
    with LazyLogging {

  // Cached catalog. Refreshed lazily when older than maxCatalogLifetime.
  private var costCatalog: ExpiringCostCatalog =
    ExpiringCostCatalog(CostCatalogService.processCostCatalog(fetchNewCatalog), java.time.LocalDateTime.now())

  /**
   * Returns the SKU for a given key, if it exists
   */
  def getSku(key: CostCatalogKey): Option[CostCatalogValue] = getOrFetchCachedCatalog().get(key)

  def serviceRegistryActor: ActorRef = serviceRegistry
  override def receive: Receive = {
    case ShutdownCommand => context stop self
    case other =>
      logger.error(
        s"Programmer Error: Unexpected message ${other.toPrettyElidedString(1000)} received by ${this.self.path.name}."
      )
  }

  private def fetchNewCatalog: Iterable[Sku] =
    mockTestData.getOrElse(CostCatalogService.makeInitialWebRequest().iterateAll().asScala)

  def getCatalogAge: Duration = Duration.between(costCatalog.fetchTime, java.time.LocalDateTime.now())
  private def isCurrentCatalogExpired: Boolean = getCatalogAge.toNanos > maxCatalogLifetime.toNanos

  private def getOrFetchCachedCatalog(): Map[CostCatalogKey, CostCatalogValue] = {
    if (isCurrentCatalogExpired) {
      costCatalog =
        ExpiringCostCatalog(CostCatalogService.processCostCatalog(fetchNewCatalog), java.time.LocalDateTime.now())
    }
    costCatalog.catalog
  }
}

object CostCatalogService {
  // Can be gleaned by using googleClient.listServices
  private val COMPUTE_ENGINE_SERVICE_NAME = "services/6F81-5844-456A"

  // ISO 4217 https://developers.google.com/adsense/management/appendix/currencies
  private val DEFAULT_CURRENCY_CODE = "USD"

  private lazy val googleClient = CloudCatalogClient.create

  /**
   * Organizes the response from google (a flat list of Sku objects) into a map that we can use for efficient lookups.
   * NB: This function takes an iterable so we can take advantage of the pagination provided by googleClient.listSkus.
   * Ideally, we don't want to have an entire, unprocessed, cost catalog in memory at once since it's ~20MB.
   */
  private def processCostCatalog(skus: Iterable[Sku]): Map[CostCatalogKey, CostCatalogValue] =
    // TODO: Account for key collisions (same key can be in multiple regions)
    // TODO: reduce memory footprint of returned map  (don't store entire SKU object)
    skus.foldLeft(Map.empty[CostCatalogKey, CostCatalogValue]) { case (acc, sku) =>
      acc + CostCatalogUtils.convertSkuToKeyValuePair(sku)
    }

  /**
   * Makes a web request to google that lists all SKUs. The request is paginated, so the response will only include
   * the first page. In order to get to all SKUs via iterateAll(), additional web requests will be made.
   */
  private def makeInitialWebRequest(): CloudCatalogClient.ListSkusPagedResponse = {
    val request = ListSkusRequest
      .newBuilder()
      .setParent(COMPUTE_ENGINE_SERVICE_NAME)
      .setCurrencyCode(DEFAULT_CURRENCY_CODE)
      .build()
    googleClient.listSkus(request)
  }

  /**
   * Helper function to save mock testing data.
   */
  def saveResponseToFileForTesting(filepath: String): Unit = {
    val dataToSave: List[Sku] = makeInitialWebRequest().iterateAll().asScala.toList
    val fileOutputStream = new FileOutputStream(filepath)
    val objectOutputStream = new ObjectOutputStream(fileOutputStream)
    objectOutputStream.writeObject(dataToSave)
    objectOutputStream.flush()
  }
}
