package cromwell.services.cost

import akka.actor.{Actor, ActorRef}
import com.google.`type`.Money
import com.google.cloud.billing.v1._
import com.typesafe.config.Config
import com.typesafe.scalalogging.LazyLogging
import common.util.StringUtil.EnhancedToStringable
import cromwell.services.ServiceRegistryActor.ServiceRegistryMessage
import cromwell.services.cost.GcpCostCatalogService.{COMPUTE_ENGINE_SERVICE_NAME, DEFAULT_CURRENCY_CODE}
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import java.time.{Duration, Instant}
import scala.jdk.CollectionConverters.IterableHasAsScala
import java.time.temporal.ChronoUnit.SECONDS
import scala.util.{Failure, Success, Try}

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

}

case class GcpCostLookupRequest(vmInfo: InstantiatedVmInfo, replyTo: ActorRef) extends ServiceRegistryMessage {
  override def serviceName: String = "GcpCostCatalogService"
}
case class GcpCostLookupResponse(calculatedCost: Option[BigDecimal])
case class CostCatalogValue(catalogObject: Sku)
case class ExpiringGcpCostCatalog(catalog: Map[CostCatalogKey, CostCatalogValue], fetchTime: Instant)

object GcpCostCatalogService {
  // Can be gleaned by using googleClient.listServices
  private val COMPUTE_ENGINE_SERVICE_NAME = "services/6F81-5844-456A"

  // ISO 4217 https://developers.google.com/adsense/management/appendix/currencies
  private val DEFAULT_CURRENCY_CODE = "USD"
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

  // See: https://cloud.google.com/billing/v1/how-tos/catalog-api
  private def calculateCpuPricePerHour(cpuSku: Sku, coreCount: Int): Try[BigDecimal] = {
    val pricingInfo = getMostRecentPricingInfo(cpuSku)
    val usageUnit = pricingInfo.getPricingExpression.getUsageUnit
    if (usageUnit != "h") {
      return Failure(new UnsupportedOperationException(s"Expected usage units of CPUs to be 'h'. Got ${usageUnit}"))
    }
    // Price per hour of a single core
    // NB: Ignoring "TieredRates" here (the idea that stuff gets cheaper the more you use).
    // Technically, we should write code that determines which tier(s) to use.
    // In practice, from what I've seen, CPU cores and RAM don't have more than a single tier.
    val costPerUnit: Money = pricingInfo.getPricingExpression.getTieredRates(0).getUnitPrice
    val costPerCorePerHour: BigDecimal =
      costPerUnit.getUnits + (costPerUnit.getNanos * 10e-9) // Same as above, but as a big decimal

    println(s"Calculated ${coreCount} cpu cores to cost ${costPerCorePerHour * coreCount} per hour")
    Success(costPerCorePerHour * coreCount)
  }

  private def calculateRamPricePerHour(ramSku: Sku, ramGbCount: Int): Try[BigDecimal] = {
    val pricingInfo = getMostRecentPricingInfo(ramSku)
    val usageUnit = pricingInfo.getPricingExpression.getUsageUnit
    if (usageUnit != "GiBy.h") {
      return Failure(new UnsupportedOperationException(s"Expected usage units of RAM to be 'GiBy.h'. Got ${usageUnit}"))
    }
    val costPerUnit: Money = pricingInfo.getPricingExpression.getTieredRates(0).getUnitPrice
    val costPerGbHour: BigDecimal =
      costPerUnit.getUnits + (costPerUnit.getNanos * 10e-9) // Same as above, but as a big decimal
    println(s"Calculated ${ramGbCount} GB of ram to cost ${ramGbCount * costPerGbHour} per hour")
    Success(costPerGbHour * ramGbCount)
  }

  private def getMostRecentPricingInfo(sku: Sku): PricingInfo = {
    val mostRecentPricingInfoIndex = sku.getPricingInfoCount - 1
    sku.getPricingInfo(mostRecentPricingInfoIndex)
  }

  private def calculateVmCostPerHour(instantiatedVmInfo: InstantiatedVmInfo): Try[BigDecimal] = {
    val machineType = MachineType.fromGoogleMachineTypeString(instantiatedVmInfo.machineType)
    val usageType = UsageType.fromBoolean(instantiatedVmInfo.preemptible)
    val machineCustomization = MachineCustomization.fromMachineTypeString(instantiatedVmInfo.machineType)
    val region = instantiatedVmInfo.region
    val coreCount = MachineType.extractCoreCountFromMachineTypeString(instantiatedVmInfo.machineType)
    val ramMbCount = MachineType.extractRamMbFromMachineTypeString(instantiatedVmInfo.machineType)
    val ramGbCount = ramMbCount.getOrElse(0) / 1024

    val cpuResourceType = Cpu
    val cpuKey =
      CostCatalogKey(machineType.get, usageType, machineCustomization, cpuResourceType, region) // TODO .get

    // As of Sept 2024 the cost catalog does not contain entries for custom N1 machines. If we're using N1, attempt
    // to fall back to predefined.
    lazy val n1PredefinedCpuSku = (machineType, machineCustomization) match {
      case (Some(N1), Custom) => getSku(cpuKey.copy(machineCustomization = Predefined))
      case _ => None
    }

    val cpuSku = getSku(cpuKey).orElse(n1PredefinedCpuSku)
    if (cpuSku.isEmpty) {
      println(s"Failed to find CPU Sku for ${cpuKey}")
    } else {
      println(s"Found CPU Sku ${cpuSku.get.catalogObject.getDescription} from key ${cpuKey}")
    }
    val cpuCost = cpuSku.map(sku => calculateCpuPricePerHour(sku.catalogObject, coreCount.get)) // TODO .get

    val ramResourceType = Ram
    val ramKey =
      CostCatalogKey(machineType.get, usageType, machineCustomization, ramResourceType, region) // TODO .get

    // As of Sept 2024 the cost catalog does not contain entries for custom N1 machines. If we're using N1, attempt
    // to fall back to predefined.
    lazy val n1PredefinedRamSku = (machineType, machineCustomization) match {
      case (Some(N1), Custom) => getSku(ramKey.copy(machineCustomization = Predefined))
      case _ => None
    }

    val ramSku = getSku(ramKey).orElse(n1PredefinedRamSku)
    if (ramSku.isEmpty) {
      println(s"Failed to find Ram Sku for ${ramKey}")
    } else {
      println(s"Found CPU Sku ${ramSku.get.catalogObject.getDescription} from key ${ramKey}")
    }
    val ramCost = ramSku.map(sku => calculateRamPricePerHour(sku.catalogObject, ramGbCount)) // TODO .get
    Success(cpuCost.get.get + ramCost.get.get)
  }

  def serviceRegistryActor: ActorRef = serviceRegistry
  override def receive: Receive = {
    case GcpCostLookupRequest(vmInfo, replyTo) =>
      val calculatedCost = calculateVmCostPerHour(vmInfo).toOption
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
