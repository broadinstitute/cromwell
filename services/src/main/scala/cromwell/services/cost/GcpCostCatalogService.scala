package cromwell.services.cost

import akka.actor.{Actor, ActorRef}
import cats.implicits.catsSyntaxValidatedId
import com.google.`type`.Money
import com.google.cloud.billing.v1._
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.util.StringUtil.EnhancedToStringable
import common.validation.ErrorOr._
import common.validation.ErrorOr.ErrorOr
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.services.cost.GcpCostCatalogService.{COMPUTE_ENGINE_SERVICE_NAME, DEFAULT_CURRENCY_CODE}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import java.time.{Duration, Instant}
import scala.jdk.CollectionConverters.IterableHasAsScala
import java.time.temporal.ChronoUnit.SECONDS

case class CostCatalogKey(machineType: MachineType,
                          usageType: UsageType,
                          machineCustomization: MachineCustomization,
                          resourceType: ResourceType,
                          region: String
)

object CostCatalogKey {

  // Specifically support only the SKUs that we know we can use. This is brittle and I hate it, but the more structured
  // fields available in the SKU don't give us enough information without relying on the human-readable descriptions.
  //
  // N1: We usually use custom machines but SKUs are only available for predefined; we'll fall back to these SKUs.
  // N2 and N2D: We only use custom machines.

  // Use this regex to filter down to just the SKUs we are interested in.
  // NB: This should be updated if we add new machine types or the cost catalog descriptions change
  final val expectedSku =
    (".*?N1 Predefined Instance (Core|Ram) .*|" +
      ".*?N2 Custom Instance (Core|Ram) .*|" +
      ".*?N2D AMD Custom Instance (Core|Ram) .*").r

  def apply(sku: Sku): List[CostCatalogKey] =
    for {
      _ <- expectedSku.findFirstIn(sku.getDescription).toList
      machineType <- MachineType.fromSku(sku).toList
      resourceType <- ResourceType.fromSku(sku).toList
      usageType <- UsageType.fromSku(sku).toList
      machineCustomization <- MachineCustomization.fromSku(sku).toList
      region <- sku.getServiceRegionsList.asScala.toList
    } yield CostCatalogKey(machineType, usageType, machineCustomization, resourceType, region)

  def apply(instantiatedVmInfo: InstantiatedVmInfo, resourceType: ResourceType): ErrorOr[CostCatalogKey] =
    MachineType.fromGoogleMachineTypeString(instantiatedVmInfo.machineType).map { mType =>
      CostCatalogKey(
        mType,
        UsageType.fromBoolean(instantiatedVmInfo.preemptible),
        MachineCustomization.fromMachineTypeString(instantiatedVmInfo.machineType),
        resourceType,
        instantiatedVmInfo.region
      )
    }
}

case class GcpCostLookupRequest(vmInfo: InstantiatedVmInfo, replyTo: ActorRef) extends ServiceRegistryMessage {
  override def serviceName: String = "GcpCostCatalogService"
}
case class GcpCostLookupResponse(calculatedCost: ErrorOr[BigDecimal])
case class CostCatalogValue(catalogObject: Sku)
case class ExpiringGcpCostCatalog(catalog: Map[CostCatalogKey, CostCatalogValue], fetchTime: Instant)

object GcpCostCatalogService {
  // Can be gleaned by using googleClient.listServices
  private val COMPUTE_ENGINE_SERVICE_NAME = "services/6F81-5844-456A"

  // ISO 4217 https://developers.google.com/adsense/management/appendix/currencies
  private val DEFAULT_CURRENCY_CODE = "USD"

  def getMostRecentPricingInfo(sku: Sku): PricingInfo = {
    val mostRecentPricingInfoIndex = sku.getPricingInfoCount - 1
    sku.getPricingInfo(mostRecentPricingInfoIndex)
  }

  // See: https://cloud.google.com/billing/v1/how-tos/catalog-api
  def calculateCpuPricePerHour(cpuSku: Sku, coreCount: Int): ErrorOr[BigDecimal] = {
    val pricingInfo = getMostRecentPricingInfo(cpuSku)
    val usageUnit = pricingInfo.getPricingExpression.getUsageUnit
    if (usageUnit == "h") {
      // Price per hour of a single core
      // NB: Ignoring "TieredRates" here (the idea that stuff gets cheaper the more you use).
      // Technically, we should write code that determines which tier(s) to use.
      // In practice, from what I've seen, CPU cores and RAM don't have more than a single tier.
      val costPerUnit: Money = pricingInfo.getPricingExpression.getTieredRates(0).getUnitPrice
      val costPerCorePerHour: BigDecimal =
        costPerUnit.getUnits + (costPerUnit.getNanos * 10e-9) // Same as above, but as a big decimal

      println(s"Calculated ${coreCount} cpu cores to cost ${costPerCorePerHour * coreCount} per hour")
      val result = costPerCorePerHour * coreCount
      result.validNel
    } else {
      s"Expected usage units of CPUs to be 'h'. Got ${usageUnit}".invalidNel
    }
  }

  def calculateRamPricePerHour(ramSku: Sku, ramGbCount: Int): ErrorOr[BigDecimal] = {
    val pricingInfo = getMostRecentPricingInfo(ramSku)
    val usageUnit = pricingInfo.getPricingExpression.getUsageUnit
    if (usageUnit == "GiBy.h") {
      val costPerUnit: Money = pricingInfo.getPricingExpression.getTieredRates(0).getUnitPrice
      val costPerGbHour: BigDecimal =
        costPerUnit.getUnits + (costPerUnit.getNanos * 10e-9) // Same as above, but as a big decimal

      println(s"Calculated ${ramGbCount} GB of ram to cost ${ramGbCount * costPerGbHour} per hour")
      val result = costPerGbHour * ramGbCount
      result.validNel
    } else {
      s"Expected usage units of RAM to be 'GiBy.h'. Got ${usageUnit}".invalidNel
    }
  }
}

/**
 * This actor handles the fetching and caching of the public Google Cost Catalog.
 * This enables us to look up the SKUs necessary for cost calculations.
 */
class GcpCostCatalogService(serviceConfig: Config, globalConfig: Config, serviceRegistry: ActorRef)
    extends Actor
    with LazyLogging {

  private val maxCatalogLifetime: Duration =
    Duration.of(CostCatalogConfig(serviceConfig).catalogExpirySeconds.longValue, SECONDS)

  private var googleClient: Option[CloudCatalogClient] = Option.empty

  // Cached catalog. Refreshed lazily when older than maxCatalogLifetime.
  private var costCatalog: Option[ExpiringGcpCostCatalog] = Option.empty

  /**
   * Returns the SKU for a given key, if it exists
   */
  def getSku(key: CostCatalogKey): Option[CostCatalogValue] = getOrFetchCachedCatalog().get(key)

  protected def fetchNewCatalog: Iterable[Sku] = {
    if (googleClient.isEmpty) {
      // We use option rather than lazy here so that the client isn't created when it is told to shutdown (see receive override)
      googleClient = Some(CloudCatalogClient.create)
    }
    makeInitialWebRequest(googleClient.get).iterateAll().asScala
  }

  def getCatalogAge: Duration =
    Duration.between(costCatalog.map(c => c.fetchTime).getOrElse(Instant.ofEpochMilli(0)), Instant.now())
  private def isCurrentCatalogExpired: Boolean = getCatalogAge.toNanos > maxCatalogLifetime.toNanos

  private def getOrFetchCachedCatalog(): Map[CostCatalogKey, CostCatalogValue] = {
    if (costCatalog.isEmpty || isCurrentCatalogExpired) {
      logger.info("Fetching a new GCP public cost catalog.")
      costCatalog = Some(ExpiringGcpCostCatalog(processCostCatalog(fetchNewCatalog), Instant.now()))
    }
    costCatalog.map(expiringCatalog => expiringCatalog.catalog).getOrElse(Map.empty)
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

  /**
   * Organizes the response from google (a flat list of Sku objects) into a map that we can use for efficient lookups.
   * NB: This function takes an iterable so we can take advantage of the pagination provided by googleClient.listSkus.
   * Ideally, we don't want to have an entire, unprocessed, cost catalog in memory at once since it's ~20MB.
   */
  private def processCostCatalog(skus: Iterable[Sku]): Map[CostCatalogKey, CostCatalogValue] =
    skus.foldLeft(Map.empty[CostCatalogKey, CostCatalogValue]) { case (acc, sku) =>
      val keys = CostCatalogKey(sku)

      // We expect that every cost catalog key is unique, but changes to the SKUs returned by Google may
      // break this assumption. Check and log an error if we find collisions.
      val collisions = keys.filter(acc.contains)
      if (collisions.nonEmpty)
        logger.error(s"Found SKU key collision for ${sku.getDescription}")

      acc ++ keys.map(k => (k, CostCatalogValue(sku)))
    }

  def lookUpSku(instantiatedVmInfo: InstantiatedVmInfo, resourceType: ResourceType): ErrorOr[Sku] =
    CostCatalogKey(instantiatedVmInfo, resourceType).flatMap { key =>
      // As of Sept 2024 the cost catalog does not contain entries for custom N1 machines. If we're using N1, attempt
      // to fall back to predefined.
      lazy val n1PredefinedKey =
        (key.machineType, key.machineCustomization) match {
          case (N1, Custom) => Option(key.copy(machineCustomization = Predefined))
          case _ => None
        }
      val sku = getSku(key).orElse(n1PredefinedKey.flatMap(getSku)).map(_.catalogObject)
      sku match {
        case Some(sku) => sku.validNel
        case None => s"Failed to look up ${resourceType} SKU for ${instantiatedVmInfo}".invalidNel
      }
    }

  def calculateVmCostPerHour(instantiatedVmInfo: InstantiatedVmInfo): ErrorOr[BigDecimal] = {
    val cpuPricePerHour: ErrorOr[BigDecimal] = for {
      cpuSku <- lookUpSku(instantiatedVmInfo, Cpu)
      coreCount <- MachineType.extractCoreCountFromMachineTypeString(instantiatedVmInfo.machineType)
      pricePerHour <- GcpCostCatalogService.calculateCpuPricePerHour(cpuSku, coreCount)
    } yield pricePerHour

    val ramPricePerHour: ErrorOr[BigDecimal] = for {
      ramSku <- lookUpSku(instantiatedVmInfo, Ram)
      ramMbCount <- MachineType.extractRamMbFromMachineTypeString(instantiatedVmInfo.machineType)
      ramGbCount = ramMbCount / 1024
      pricePerHour <- GcpCostCatalogService.calculateRamPricePerHour(ramSku, ramGbCount)
    } yield pricePerHour

    for {
      cpuPrice <- cpuPricePerHour
      ramPrice <- ramPricePerHour
    } yield cpuPrice + ramPrice
  }

  def serviceRegistryActor: ActorRef = serviceRegistry
  override def receive: Receive = {
    case GcpCostLookupRequest(vmInfo, replyTo) =>
      val calculatedCost = calculateVmCostPerHour(vmInfo)
      val response = GcpCostLookupResponse(calculatedCost)
      replyTo ! response
    case ShutdownCommand =>
      googleClient.foreach(client => client.shutdownNow())
      context stop self
    case other =>
      logger.error(
        s"Programmer Error: Unexpected message ${other.toPrettyElidedString(1000)} received by ${this.self.path.name}."
      )
  }
}
