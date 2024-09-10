package cromwell.backend.google.pipelines.common

import akka.testkit.{ImplicitSender, TestActorRef, TestProbe}
import cromwell.core.TestKitSuite
import org.scalatest.flatspec.{AnyFlatSpecLike}
import org.scalatest.matchers.should.Matchers
import cromwell.services.cost._
import org.scalatest.concurrent.Eventually

class CostLookupSpec
  extends TestKitSuite
    with AnyFlatSpecLike
    with Matchers
    with Eventually
    with ImplicitSender {
  behavior of "CostLookup"

  def constructTestActor: GcpCostCatalogServiceTestActor =
    TestActorRef(
      new GcpCostCatalogServiceTestActor(GcpCostCatalogServiceSpec.config,
        GcpCostCatalogServiceSpec.config,
        TestProbe().ref
      )
    ).underlyingActor

  val testCatalogService = constructTestActor

  it should "find a CPU from Runtime Attributes" in {
    val machineType = Some(N2)
    val usageType = Some(OnDemand)
    val customization = Some(Custom)
    val resourceGroup = Some(Cpu)
    val region = "europe-west9"
    val key = CostCatalogKey(machineType, usageType, customization, resourceGroup, region)
    val result = testCatalogService.getSku(key)
    result shouldBe true
  }
}
