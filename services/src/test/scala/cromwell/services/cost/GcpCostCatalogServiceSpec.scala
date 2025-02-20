package cromwell.services.cost

import akka.actor.ActorRef
import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import com.google.cloud.billing.v1.{CloudCatalogClient, Sku}
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import org.scalatest.prop.TableDrivenPropertyChecks

import java.time.Duration
import java.io.{File, FileInputStream, FileOutputStream}
import scala.collection.mutable.ListBuffer
import scala.util.Using

object GcpCostCatalogServiceSpec {
  val catalogExpirySeconds: Long = 1 // Short duration so we can do a cache expiry test
  val config: Config = ConfigFactory.parseString(
    s"catalogExpirySeconds = $catalogExpirySeconds, enabled = true"
  )
  val mockTestDataFilePath: String = "services/src/test/scala/cromwell/services/cost/serializedSkuList.testData"
}
class GcpCostCatalogServiceTestActor(serviceConfig: Config, globalConfig: Config, serviceRegistry: ActorRef)
    extends GcpCostCatalogService(serviceConfig, globalConfig, serviceRegistry) {
  def loadMockData: Iterable[Sku] = {
    val fis = new FileInputStream(new File(GcpCostCatalogServiceSpec.mockTestDataFilePath))
    val skuList = ListBuffer[Sku]()
    try {
      var sku: Sku = null
      do {
        sku = Sku.parseDelimitedFrom(fis)
        if (sku != null) skuList += sku
      } while (sku != null)
    } finally
      fis.close()
    skuList.toSeq
  }
  def saveMockData(): Unit =
    Using.resources(CloudCatalogClient.create,
                    new FileOutputStream(new File(GcpCostCatalogServiceSpec.mockTestDataFilePath))
    ) { (googleClient, fos) =>
      val skus = super.fetchSkuIterable(googleClient)
      skus.foreach { sku =>
        sku.writeDelimitedTo(fos)
      }
    }

  override def receive: Receive = { case ShutdownCommand =>
    context stop self
  }

  override def fetchNewCatalog: ExpiringGcpCostCatalog = makeCatalog(loadMockData)
}

class GcpCostCatalogServiceSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with Eventually
    with ImplicitSender
    with TableDrivenPropertyChecks {
  behavior of "CostCatalogService"

  def constructTestActor: GcpCostCatalogServiceTestActor =
    TestActorRef(
      new GcpCostCatalogServiceTestActor(GcpCostCatalogServiceSpec.config,
                                         GcpCostCatalogServiceSpec.config,
                                         TestProbe().ref
      )
    ).underlyingActor
  private val testActorRef = constructTestActor
  private val minimumExpectedSkus = 10

  it should "properly load mock test data file" in {
    testActorRef.loadMockData.toList.length shouldBe >(minimumExpectedSkus)
  }

  it should "cache catalogs properly" in {
    val testLookupKey = CostCatalogKey(
      resourceInfo = N2,
      usageType = Preemptible,
      machineCustomization = Some(Predefined),
      resourceType = Cpu,
      region = "europe-west9"
    )

    val freshActor = constructTestActor
    val shortDuration: Duration = Duration.ofSeconds(GcpCostCatalogServiceSpec.catalogExpirySeconds)

    // Sanity check that the catalog has been initially fetched with a new timestamp
    freshActor.getSku(testLookupKey) // do a lookup to trigger any lazy loading
    freshActor.getCatalogAge.toNanos should (be > 0.toLong)
    freshActor.getCatalogAge.toNanos should (be < shortDuration.toNanos)

    // Simulate the cached catalog living longer than its lifetime
    Thread.sleep(shortDuration.plus(shortDuration).toMillis)

    // Confirm that the catalog is old
    freshActor.getCatalogAge.toNanos should (be > shortDuration.toNanos)

    // Check that doing a lookup causes the catalog to be refreshed
    freshActor.getSku(testLookupKey)
    freshActor.getCatalogAge.toNanos should (be > 0.toLong)
    freshActor.getCatalogAge.toNanos should (be < shortDuration.toNanos)
  }

  it should "find CPU and RAM skus for all supported machine types" in {
    val lookupRows = Table(
      ("machineType", "usage", "customization", "resource", "region", "exists"),
      (N1, Preemptible, Some(Predefined), Cpu, "us-west1", true),
      (N1, Preemptible, Some(Predefined), Ram, "us-west1", true),
      (N1, OnDemand, Some(Predefined), Cpu, "us-west1", true),
      (N1, OnDemand, Some(Predefined), Ram, "us-west1", true),
      (N1, Preemptible, Some(Custom), Cpu, "us-west1", false),
      (N1, Preemptible, Some(Custom), Ram, "us-west1", false),
      (N1, OnDemand, Some(Custom), Cpu, "us-west1", false),
      (N1, OnDemand, Some(Custom), Ram, "us-west1", false),
      (N2, Preemptible, Some(Predefined), Cpu, "us-west1", false),
      (N2, Preemptible, Some(Predefined), Ram, "us-west1", false),
      (N2, OnDemand, Some(Predefined), Cpu, "us-west1", false),
      (N2, OnDemand, Some(Predefined), Ram, "us-west1", false),
      (N2, Preemptible, Some(Custom), Cpu, "us-west1", true),
      (N2, Preemptible, Some(Custom), Ram, "us-west1", true),
      (N2, OnDemand, Some(Custom), Cpu, "us-west1", true),
      (N2, OnDemand, Some(Custom), Ram, "us-west1", true),
      (N2d, Preemptible, Some(Predefined), Cpu, "us-west1", false),
      (N2d, Preemptible, Some(Predefined), Ram, "us-west1", false),
      (N2d, OnDemand, Some(Predefined), Cpu, "us-west1", false),
      (N2d, OnDemand, Some(Predefined), Ram, "us-west1", false),
      (N2d, Preemptible, Some(Custom), Cpu, "us-west1", true),
      (N2d, Preemptible, Some(Custom), Ram, "us-west1", true),
      (N2d, OnDemand, Some(Custom), Cpu, "us-west1", true),
      (N2d, OnDemand, Some(Custom), Ram, "us-west1", true)
    )

    forAll(lookupRows) { case (machineType, usage, customization, resource, region, exists: Boolean) =>
      val key = CostCatalogKey(machineType, usage, customization, resource, region)
      val result = testActorRef.getSku(key)
      result.nonEmpty shouldBe exists
    }
  }

  it should "find the skus for a VM when appropriate" in {
    val lookupRows = Table(
      ("instantiatedVmInfo", "resource", "skuDescription"),
      (InstantiatedVmInfo("europe-west9", "custom-16-32768", None, false),
       Cpu,
       "N1 Predefined Instance Core running in Paris"
      ),
      (InstantiatedVmInfo("europe-west9", "custom-16-32768", None, false),
       Ram,
       "N1 Predefined Instance Ram running in Paris"
      ),
      (InstantiatedVmInfo("us-central1", "custom-4-4096", None, true),
       Cpu,
       "Spot Preemptible N1 Predefined Instance Core running in Americas"
      ),
      (InstantiatedVmInfo("us-central1", "custom-4-4096", None, true),
       Ram,
       "Spot Preemptible N1 Predefined Instance Ram running in Americas"
      ),
      (InstantiatedVmInfo("europe-west9", "n1-custom-16-32768", None, false),
       Cpu,
       "N1 Predefined Instance Core running in Paris"
      ),
      (InstantiatedVmInfo("europe-west9", "n1-custom-16-32768", None, false),
       Ram,
       "N1 Predefined Instance Ram running in Paris"
      ),
      (InstantiatedVmInfo("us-central1", "n1-custom-4-4096", None, true),
       Cpu,
       "Spot Preemptible N1 Predefined Instance Core running in Americas"
      ),
      (InstantiatedVmInfo("us-central1", "n1-custom-4-4096", None, true),
       Ram,
       "Spot Preemptible N1 Predefined Instance Ram running in Americas"
      ),
      (InstantiatedVmInfo("us-central1", "n2-custom-4-4096", None, true),
       Cpu,
       "Spot Preemptible N2 Custom Instance Core running in Americas"
      ),
      (InstantiatedVmInfo("us-central1", "n2-custom-4-4096", None, true),
       Ram,
       "Spot Preemptible N2 Custom Instance Ram running in Americas"
      ),
      (InstantiatedVmInfo("us-central1", "n2-custom-4-4096", None, false),
       Cpu,
       "N2 Custom Instance Core running in Americas"
      ),
      (InstantiatedVmInfo("us-central1", "n2-custom-4-4096", None, false),
       Ram,
       "N2 Custom Instance Ram running in Americas"
      ),
      (InstantiatedVmInfo("us-central1", "n2d-custom-4-4096", None, true),
       Cpu,
       "Spot Preemptible N2D AMD Custom Instance Core running in Americas"
      ),
      (InstantiatedVmInfo("us-central1", "n2d-custom-4-4096", None, true),
       Ram,
       "Spot Preemptible N2D AMD Custom Instance Ram running in Americas"
      ),
      (InstantiatedVmInfo("us-central1", "n2d-custom-4-4096", None, false),
       Cpu,
       "N2D AMD Custom Instance Core running in Americas"
      ),
      (InstantiatedVmInfo("us-central1", "n2d-custom-4-4096", None, false),
       Ram,
       "N2D AMD Custom Instance Ram running in Americas"
      )
    )

    forAll(lookupRows) { case (instantiatedVmInfo: InstantiatedVmInfo, resource: ResourceType, expectedSku: String) =>
      val skuOr = testActorRef.lookUpSku(instantiatedVmInfo, resource)
      skuOr.isValid shouldBe true
      skuOr.map(sku => sku.getDescription shouldEqual expectedSku)
    }
  }

  it should "fail to find the skus for a VM when appropriate" in {
    val lookupRows = Table(
      ("instantiatedVmInfo", "resource", "errors"),
      (InstantiatedVmInfo("us-central1", "custooooooom-4-4096", None, true),
       Cpu,
       List("Unrecognized machine type: custooooooom-4-4096")
      ),
      (InstantiatedVmInfo("us-central1", "n2custom-4-4096", None, true),
       Cpu,
       List("Unrecognized machine type: n2custom-4-4096")
      ),
      (InstantiatedVmInfo("us-central1", "standard-4-4096", None, true),
       Cpu,
       List("Unrecognized machine type: standard-4-4096")
      ),
      (InstantiatedVmInfo("planet-mars1", "custom-4-4096", None, true),
       Cpu,
       List("Failed to look up Cpu SKU for InstantiatedVmInfo(planet-mars1,custom-4-4096,None,true)")
      )
    )

    forAll(lookupRows) {
      case (instantiatedVmInfo: InstantiatedVmInfo, resource: ResourceType, expectedErrors: List[String]) =>
        val skuOr = testActorRef.lookUpSku(instantiatedVmInfo, resource)
        skuOr.isValid shouldBe false
        skuOr.leftMap(errors => errors.toList shouldEqual expectedErrors)
    }
  }

  it should "calculate the cost per hour for a VM" in {
    // Create BigDecimals from strings to avoid inequality due to floating point shenanigans
    val lookupRows = Table(
      ("instantiatedVmInfo", "costPerHour"),
      (InstantiatedVmInfo("us-central1", "custom-4-4096", None, true), BigDecimal(".0361")),
      (InstantiatedVmInfo("us-central1", "n2-custom-4-4096", None, true), BigDecimal(".04254400000000000480")),
      (InstantiatedVmInfo("us-central1", "n2d-custom-4-4096", None, true), BigDecimal(".02371600000000000040")),
      (InstantiatedVmInfo("us-central1", "custom-4-4096", None, false), BigDecimal(".143392")),
      (InstantiatedVmInfo("us-central1", "n2-custom-4-4096", None, false), BigDecimal(".150561600")),
      (InstantiatedVmInfo("us-central1", "n2d-custom-4-4096", None, false), BigDecimal(".130989600000000012")),
      (InstantiatedVmInfo("europe-west9", "custom-4-4096", None, true), BigDecimal(".035018080000000004")),
      (InstantiatedVmInfo("europe-west9", "n2-custom-4-4096", None, true), BigDecimal("0.049532000000000004")),
      (InstantiatedVmInfo("europe-west9", "n2d-custom-4-4096", None, true), BigDecimal("0.030608000000000004")),
      (InstantiatedVmInfo("europe-west9", "custom-4-4096", None, false), BigDecimal(".1663347200000000040")),
      (InstantiatedVmInfo("europe-west9", "n2-custom-4-4352", None, false), BigDecimal(".17594163050")),
      (InstantiatedVmInfo("europe-west9", "n2d-custom-4-4096", None, false), BigDecimal(".151947952"))
    )

    forAll(lookupRows) { case (instantiatedVmInfo: InstantiatedVmInfo, expectedCostPerHour: BigDecimal) =>
      val costOr = testActorRef.calculateVmCostPerHour(instantiatedVmInfo)
      costOr.isValid shouldBe true
      costOr.map(cost => cost shouldEqual expectedCostPerHour)
    }
  }

  it should "fail to calculate the cost oer hour for a VM" in {

    val lookupRows = Table(
      ("instantiatedVmInfo", "errors"),
      (InstantiatedVmInfo("us-central1", "custooooooom-4-4096", None, true),
       List("Unrecognized machine type: custooooooom-4-4096")
      ),
      (InstantiatedVmInfo("us-central1", "n2_custom_4_4096", None, true),
       List("Unrecognized machine type: n2_custom_4_4096")
      ),
      (InstantiatedVmInfo("us-central1", "custom-foo-4096", None, true),
       List("Could not extract core count from custom-foo-4096")
      ),
      (InstantiatedVmInfo("us-central1", "custom-16-bar", None, true),
       List("Could not extract Ram MB count from custom-16-bar")
      ),
      (InstantiatedVmInfo("us-central1", "123-456-789", None, true), List("Unrecognized machine type: 123-456-789")),
      (InstantiatedVmInfo("us-central1", "n2-16-4096", None, true),
       List("Failed to look up Cpu SKU for InstantiatedVmInfo(us-central1,n2-16-4096,None,true)")
      ),
      (InstantiatedVmInfo("planet-mars1", "n2-custom-4-4096", None, true),
       List("Failed to look up Cpu SKU for InstantiatedVmInfo(planet-mars1,n2-custom-4-4096,None,true)")
      )
    )

    forAll(lookupRows) { case (instantiatedVmInfo: InstantiatedVmInfo, expectedErrors: List[String]) =>
      val costOr = testActorRef.calculateVmCostPerHour(instantiatedVmInfo)
      costOr.isValid shouldBe false
      costOr.leftMap(errors => errors.toList shouldEqual expectedErrors)
    }
  }
}
