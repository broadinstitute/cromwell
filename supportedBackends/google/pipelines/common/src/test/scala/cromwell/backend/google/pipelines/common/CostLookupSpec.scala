package cromwell.backend.google.pipelines.common

import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import cromwell.core.TestKitSuite
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import cromwell.services.cost._
import org.scalatest.concurrent.Eventually
import org.scalatest.prop.TableDrivenPropertyChecks

class CostLookupSpec
    extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with Eventually
    with ImplicitSender
    with TableDrivenPropertyChecks {
  behavior of "CostLookup"

  def constructTestActor: GcpCostCatalogServiceTestActor =
    TestActorRef(
      new GcpCostCatalogServiceTestActor(GcpCostCatalogServiceSpec.config,
                                         GcpCostCatalogServiceSpec.config,
                                         TestProbe().ref
      )
    ).underlyingActor

  val testCatalogService = constructTestActor

  it should "find a CPU sku" in {
    val machineType = N2
    val usageType = OnDemand
    val customization = Custom
    val resourceGroup = Cpu
    val region = "europe-west9"
    val key = CostCatalogKey(machineType, usageType, customization, resourceGroup, region)
    val result = testCatalogService.getSku(key).get.catalogObject.getDescription
    result shouldBe "N2 Custom Instance Core running in Paris"
  }

  it should "find a RAM sku" in {
    val machineType = N2
    val usageType = OnDemand
    val customization = Custom
    val resourceGroup = Ram
    val region = "europe-west9"
    val key = CostCatalogKey(machineType, usageType, customization, resourceGroup, region)
    val result = testCatalogService.getSku(key).get.catalogObject.getDescription
    result shouldBe "N2 Custom Instance Ram running in Paris"
  }

  it should "find CPU skus for all supported machine types" in {
    val lookupRows = Table(
      ("machineType", "usage", "customization", "resource", "region", "exists"),
      (N1, Preemptible, Predefined, Cpu, "us-west1", true),
      (N1, Preemptible, Predefined, Ram, "us-west1", true),
      (N1, OnDemand, Predefined, Cpu, "us-west1", true),
      (N1, OnDemand, Predefined, Ram, "us-west1", true),
      (N1, Preemptible, Custom, Cpu, "us-west1", false),
      (N1, Preemptible, Custom, Ram, "us-west1", false),
      (N1, OnDemand, Custom, Cpu, "us-west1", false),
      (N1, OnDemand, Custom, Ram, "us-west1", false),
      (N2, Preemptible, Predefined, Cpu, "us-west1", false),
      (N2, Preemptible, Predefined, Ram, "us-west1", false),
      (N2, OnDemand, Predefined, Cpu, "us-west1", false),
      (N2, OnDemand, Predefined, Ram, "us-west1", false),
      (N2, Preemptible, Custom, Cpu, "us-west1", true),
      (N2, Preemptible, Custom, Ram, "us-west1", true),
      (N2, OnDemand, Custom, Cpu, "us-west1", true),
      (N2, OnDemand, Custom, Ram, "us-west1", true),
      (N2d, Preemptible, Predefined, Cpu, "us-west1", false),
      (N2d, Preemptible, Predefined, Ram, "us-west1", false),
      (N2d, OnDemand, Predefined, Cpu, "us-west1", false),
      (N2d, OnDemand, Predefined, Ram, "us-west1", false),
      (N2d, Preemptible, Custom, Cpu, "us-west1", true),
      (N2d, Preemptible, Custom, Ram, "us-west1", true),
      (N2d, OnDemand, Custom, Cpu, "us-west1", true),
      (N2d, OnDemand, Custom, Ram, "us-west1", true)
    )

    forAll(lookupRows) { case (machineType, usage, customization, resource, region, exists: Boolean) =>
      val key = CostCatalogKey(machineType, usage, customization, resource, region)
      val result = testCatalogService.getSku(key)
      result.nonEmpty shouldBe exists
    }
  }
}
