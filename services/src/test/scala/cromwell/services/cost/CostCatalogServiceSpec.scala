package cromwell.services.cost

import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import com.google.cloud.billing.v1.Sku
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

import java.time.Duration
import java.io.{FileInputStream, ObjectInputStream}


object CostCatalogServiceSpec {
  val Config = ConfigFactory.parseString("control-frequency = 1 second")
}

class CostCatalogServiceSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with Eventually
    with ImplicitSender {
  behavior of "CostCatalogService"

  val mockTestDataFilePath = "services/src/test/scala/cromwell/services/cost/serializedSkuList.testData"
  val minimumExpectedSkus = 10
  lazy val testActorRef = TestActorRef(
    new CostCatalogService(CostCatalogServiceSpec.Config,
                           CostCatalogServiceSpec.Config,
                           TestProbe().ref,
                           Some(loadMockTestDataFromFile())
    )
  ).underlyingActor

  /*
   * This test will download a real cost catalog and save it to disk so we can use it for mocking.
   * Only necessary to rerun if Google changes the catalog structure.
  it should "download fresh test data" in {
    val costCatalogActor = TestActorRef(new CostCatalogService(CostCatalogServiceSpec.Config, CostCatalogServiceSpec.Config, TestProbe().ref))
    costCatalogActor.underlyingActor.saveResponseToFileForTesting(mockTestDataFilePath)
    val sanityCheckFile = new File(mockTestDataFilePath)
    sanityCheckFile.exists() shouldBe true
  }
   */

  def loadMockTestDataFromFile(): List[Sku] = {
    val fis = new FileInputStream(mockTestDataFilePath)
    val ois = new ObjectInputStream(fis)
    ois.readObject().asInstanceOf[List[Sku]]
  }

  it should "properly load mock test data" in {
    loadMockTestDataFromFile().length shouldBe >(minimumExpectedSkus)
  }

  it should "cache catalogs properly" in {
    val shortDuration = Duration.ofMillis(10)
    val mockActorWithCustomExpiration = TestActorRef(
      new CostCatalogService(CostCatalogServiceSpec.Config,
                             CostCatalogServiceSpec.Config,
                             TestProbe().ref,
                             Some(loadMockTestDataFromFile()),
                             shortDuration
      )
    ).underlyingActor
    val testLookupKey = CostCatalogKey(
      machineType = Some(N2),
      usageType = Some(Preemptible),
      machineCustomization = Some(Predefined),
      resourceGroup = Some(Cpu)
    )

    // Sanity check that the catalog has been initially fetched with a new timestamp
    mockActorWithCustomExpiration.getSku(testLookupKey) // do a lookup to trigger any lazy loading
    mockActorWithCustomExpiration.getCatalogAge.getNano shouldBe >(0)
    mockActorWithCustomExpiration.getCatalogAge.getNano shouldBe <(shortDuration.getNano)

    // Simulate the cached catalog living past its lifetime
    Thread.sleep(shortDuration.toMillis)

    // Check that the catalog is too old
    mockActorWithCustomExpiration.getCatalogAge.getNano shouldBe >(shortDuration.getNano)

    // Check that doing a lookup causes the catalog to be refreshed
    mockActorWithCustomExpiration.getSku(testLookupKey)
    mockActorWithCustomExpiration.getCatalogAge.getNano shouldBe >(0)
    mockActorWithCustomExpiration.getCatalogAge.getNano shouldBe <(shortDuration.getNano)
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
