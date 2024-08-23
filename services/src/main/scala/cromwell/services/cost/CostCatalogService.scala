package cromwell.services.cost

import akka.actor.{Actor, ActorRef}
import com.google.cloud.billing.v1.{CloudCatalogClient, _}
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.util.StringUtil.EnhancedToStringable
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import java.time.{Duration, Instant}
import scala.jdk.CollectionConverters.IterableHasAsScala
import java.time.temporal.ChronoUnit.SECONDS

case class ExpiringGcpCostCatalog(catalog: Map[CostCatalogKey, CostCatalogValue], fetchTime: Instant)

/**
 * This actor handles the fetching and caching of the public Google Cost Catalog.
 * This enables us to look up the SKUs necessary for cost calculations.
 */
class CostCatalogService(serviceConfig: Config, globalConfig: Config, serviceRegistry: ActorRef)
    extends Actor
    with LazyLogging {

  private val maxCatalogLifetime: Duration =
    Duration.of(CostCatalogConfig(serviceConfig).catalogExpirySeconds.longValue, SECONDS)

  private lazy val googleClient: CloudCatalogClient = CloudCatalogClient.create

  // Cached catalog. Refreshed lazily when older than maxCatalogLifetime.
  private var costCatalog: Option[ExpiringGcpCostCatalog] = Option.empty

  /**
   * Returns the SKU for a given key, if it exists
   */
  def getSku(key: CostCatalogKey): Option[CostCatalogValue] = getOrFetchCachedCatalog().get(key)

  def serviceRegistryActor: ActorRef = serviceRegistry
  override def receive: Receive = {
    case ShutdownCommand =>
      googleClient.shutdownNow()
      context stop self
    case other =>
      logger.error(
        s"Programmer Error: Unexpected message ${other.toPrettyElidedString(1000)} received by ${this.self.path.name}."
      )
  }

  protected def fetchNewCatalog: Iterable[Sku] =
    CostCatalogService.makeInitialWebRequest(googleClient).iterateAll().asScala
  def getCatalogAge: Duration =
    Duration.between(costCatalog.map(c => c.fetchTime).getOrElse(Instant.ofEpochMilli(0)), Instant.now())
  private def isCurrentCatalogExpired: Boolean = getCatalogAge.toNanos > maxCatalogLifetime.toNanos

  private def getOrFetchCachedCatalog(): Map[CostCatalogKey, CostCatalogValue] = {
    if (costCatalog.isEmpty || isCurrentCatalogExpired) {
      logger.info("Fetching a new GCP public cost catalog.")
      costCatalog = Some(ExpiringGcpCostCatalog(CostCatalogService.processCostCatalog(fetchNewCatalog), Instant.now()))
    }
    costCatalog.map(expiringCatalog => expiringCatalog.catalog).getOrElse(Map.empty)
  }
}

object CostCatalogService {
  // Can be gleaned by using googleClient.listServices
  private val COMPUTE_ENGINE_SERVICE_NAME = "services/6F81-5844-456A"

  // ISO 4217 https://developers.google.com/adsense/management/appendix/currencies
  private val DEFAULT_CURRENCY_CODE = "USD"

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
  private def makeInitialWebRequest(googleClient: CloudCatalogClient): CloudCatalogClient.ListSkusPagedResponse = {
    val request = ListSkusRequest
      .newBuilder()
      .setParent(COMPUTE_ENGINE_SERVICE_NAME)
      .setCurrencyCode(DEFAULT_CURRENCY_CODE)
      .build()
    googleClient.listSkus(request)
  }
}
