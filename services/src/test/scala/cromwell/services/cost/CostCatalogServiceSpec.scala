package cromwell.services.cost

import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import com.google.cloud.billing.v1.Sku
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

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

  it should "populate test data" in {
    val costCatalogActor = TestActorRef(new CostCatalogService(CostCatalogServiceSpec.Config, CostCatalogServiceSpec.Config, TestProbe().ref))
    costCatalogActor.underlyingActor.saveResponseToFileForTesting(mockTestDataFilePath)
  }

  it should "properly load mock test data" in {
    val fis = new FileInputStream(mockTestDataFilePath)
    val ois = new ObjectInputStream(fis)
    val loadedFile: List[Sku] =  ois.readObject().asInstanceOf[List[Sku]]
    loadedFile.foreach(println)
  }

  it should "boot up" in {
    val costCatalogActor = TestActorRef(new CostCatalogService(CostCatalogServiceSpec.Config, CostCatalogServiceSpec.Config, TestProbe().ref))
    costCatalogActor.underlyingActor.getCatalog.isEmpty shouldBe false
  }
}