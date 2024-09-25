package cromwell.backend.google.pipelines.common

import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import cromwell.core.TestKitSuite
import org.scalatest.flatspec.AnyFlatSpecLike
import org.scalatest.matchers.should.Matchers
import cromwell.services.cost._
import org.scalatest.concurrent.Eventually

class CostLookupSpec extends TestKitSuite with AnyFlatSpecLike with Matchers with Eventually with ImplicitSender {
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
    val machineType = Some(N2)
    val usageType = Some(OnDemand)
    val customization = Some(Custom)
    val resourceGroup = Some(Cpu)
    val region = "europe-west9"
    val key = CostCatalogKey(machineType, usageType, customization, resourceGroup, region)
    val result = testCatalogService.getSku(key).get.catalogObject.getDescription
    result shouldBe "N2 Custom Instance Core running in Paris"
  }

  it should "find a RAM sku" in {
    val machineType = Some(N2)
    val usageType = Some(OnDemand)
    val customization = Some(Custom)
    val resourceGroup = Some(Ram)
    val region = "europe-west9"
    val key = CostCatalogKey(machineType, usageType, customization, resourceGroup, region)
    val result = testCatalogService.getSku(key).get.catalogObject.getDescription
    result shouldBe "N2 Custom Instance Ram running in Paris"
  }

  it should "find CPU skus for all supported machine types" in {
    val legalMachineTypes: List[MachineType] = List(N1, N2, N2d)
    val legalUsageTypes: List[UsageType] = List(Preemptible, OnDemand)
    val legalCustomizations: List[MachineCustomization] = List(Custom, Predefined)
    val resourceGroup: Option[ResourceType] = Some(Cpu)
    val region = "us-west1"
    for (machineType <- legalMachineTypes)
      for (usageType <- legalUsageTypes)
        for (customization <- legalCustomizations) {
          val key = CostCatalogKey(Some(machineType), Some(usageType), Some(customization), resourceGroup, region)
          val result = testCatalogService.getSku(key)
          if (!result.isEmpty) {
            println("Success")
          }
          result.isEmpty shouldBe false
        }
  }
}
