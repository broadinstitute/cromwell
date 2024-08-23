package cromwell.services.cost

import akka.actor.ActorRef
import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import com.google.cloud.billing.v1.Sku
import com.typesafe.config.{Config, ConfigFactory}
import cromwell.core.TestKitSuite
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import java.time.Duration
import java.io.{FileInputStream, ObjectInputStream}

object CostCatalogServiceSpec {
  val catalogExpirySeconds: Long = 1 // Short so we can do a cache expiry test
  val config = ConfigFactory.parseString(s"catalogExpirySeconds = $catalogExpirySeconds")
  val mockTestDataFilePath = "services/src/test/scala/cromwell/services/cost/serializedSkuList.testData"
  def loadMockDataFromFile: Iterable[Sku] = {
    val fis = new FileInputStream(CostCatalogServiceSpec.mockTestDataFilePath)
    val ois = new ObjectInputStream(fis)
    ois.readObject().asInstanceOf[List[Sku]]
  }
}

class CostCatalogServiceTestActor(serviceConfig: Config, globalConfig: Config, serviceRegistry: ActorRef)
    extends CostCatalogService(serviceConfig, globalConfig, serviceRegistry) {
  override def fetchNewCatalog: Iterable[Sku] = CostCatalogServiceSpec.loadMockDataFromFile
}

class CostCatalogServiceSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with Eventually
    with ImplicitSender {
  behavior of "CostCatalogService"

  def constructTestActor: CostCatalogService =
    TestActorRef(
      new CostCatalogService(CostCatalogServiceSpec.config, CostCatalogServiceSpec.config, TestProbe().ref)
    ).underlyingActor
  private val testActorRef = constructTestActor
  private val minimumExpectedSkus = 10

  it should "properly load mock test data file" in {
    CostCatalogServiceSpec.loadMockDataFromFile.toList.length shouldBe >(minimumExpectedSkus)
  }

  it should "cache catalogs properly" in {
    val testLookupKey = CostCatalogKey(
      machineType = Some(N2),
      usageType = Some(Preemptible),
      machineCustomization = Some(Predefined),
      resourceGroup = Some(Cpu)
    )

    val freshActor = constructTestActor
    val shortDuration: Duration = Duration.ofSeconds(CostCatalogServiceSpec.catalogExpirySeconds)

    // Sanity check that the catalog has been initially fetched with a new timestamp
    freshActor.getSku(testLookupKey) // do a lookup to trigger any lazy loading
    freshActor.getCatalogAge.toNanos should (be > 0.toLong)
    freshActor.getCatalogAge.toNanos should (be < shortDuration.toNanos)

    // Simulate the cached catalog living longer than its lifetime
    Thread.sleep(shortDuration.toMillis)

    // Confirm that the catalog is old
    freshActor.getCatalogAge.toNanos should (be > shortDuration.toNanos)

    // Check that doing a lookup causes the catalog to be refreshed
    freshActor.getSku(testLookupKey)
    freshActor.getCatalogAge.toNanos should (be > 0.toLong)
    freshActor.getCatalogAge.toNanos should (be < shortDuration.toNanos)
  }

  it should "contain an expected SKU" in {
    val expectedKey = CostCatalogKey(
      machineType = Some(N2d),
      usageType = Some(Preemptible),
      machineCustomization = None,
      resourceGroup = Some(Ram)
    )

    val foundValue = testActorRef.getSku(expectedKey)
    foundValue.get.catalogObject.getDescription shouldBe "Spot Preemptible N2D AMD Instance Ram running in Johannesburg"
  }
}
