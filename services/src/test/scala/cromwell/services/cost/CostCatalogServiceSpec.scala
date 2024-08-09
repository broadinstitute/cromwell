package cromwell.services.cost

import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import com.typesafe.config.ConfigFactory
import cromwell.core.TestKitSuite
import org.scalatest.concurrent.Eventually
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers

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

  it should "boot up" in {
    val costCatalogActor = TestActorRef(new CostCatalogService(CostCatalogServiceSpec.Config, CostCatalogServiceSpec.Config, TestProbe().ref))
    costCatalogActor.underlyingActor.fetchPublicCostCatalog() shouldBe List()
  }
}