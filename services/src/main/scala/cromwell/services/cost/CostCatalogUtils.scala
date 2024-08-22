package cromwell.services.cost

import com.google.cloud.billing.v1.Sku

/*
 * These types reflect hardcoded strings found in a google cost catalog.
 */
sealed trait MachineType { def machineTypeName: String }
case object N1 extends MachineType { override val machineTypeName = "N1" }
case object N2 extends MachineType { override val machineTypeName = "N2" }
case object N2d extends MachineType { override val machineTypeName = "N2D" }

sealed trait UsageType { def typeName: String }
case object OnDemand extends UsageType { override val typeName = "OnDemand" }
case object Preemptible extends UsageType { override val typeName = "Preemptible" }

sealed trait MachineCustomization { def customizationName: String }
case object Custom extends MachineCustomization { override val customizationName = "Custom" }
case object Predefined extends MachineCustomization { override val customizationName = "Predefined" }

sealed trait ResourceGroup { def groupName: String }
case object Cpu extends ResourceGroup { override val groupName = "CPU" }
case object Ram extends ResourceGroup { override val groupName = "RAM" }
case object N1Standard extends ResourceGroup { override val groupName = "N1Standard" }

case class CostCatalogKey(machineType: Option[MachineType],
                          usageType: Option[UsageType],
                          machineCustomization: Option[MachineCustomization],
                          resourceGroup: Option[ResourceGroup]
)

case class CostCatalogValue(catalogObject: Sku)

/**
 * Utils for converting google Sku objects into a smaller more searchable format
 */
object CostCatalogUtils {
  def convertSkuToKeyValuePair(sku: Sku): (CostCatalogKey, CostCatalogValue) =
    CostCatalogKey(
      machineType = extractMachineType(sku),
      usageType = extractUsageType(sku),
      machineCustomization = extractMachineCustomization(sku),
      resourceGroup = extractResourceGroup(sku)
    ) -> CostCatalogValue(sku)

  private def extractMachineType(sku: Sku): Option[MachineType] = {
    val tokenizedDescription = sku.getDescription.split(" ")
    if (tokenizedDescription.contains(N1.machineTypeName)) Some(N1)
    else if (tokenizedDescription.contains(N2.machineTypeName)) Some(N2)
    else if (tokenizedDescription.contains(N2d.machineTypeName)) Some(N2d)
    else Option.empty
  }

  private def extractUsageType(sku: Sku): Option[UsageType] =
    sku.getCategory.getUsageType match {
      case OnDemand.typeName => Some(OnDemand)
      case Preemptible.typeName => Some(Preemptible)
      case _ => Option.empty
    }

  private def extractMachineCustomization(sku: Sku): Option[MachineCustomization] = {
    val tokenizedDescription = sku.getDescription.split(" ")
    if (tokenizedDescription.contains(Predefined.customizationName)) Some(Predefined)
    else if (tokenizedDescription.contains(Custom.customizationName)) Some(Custom)
    else Option.empty
  }

  private def extractResourceGroup(sku: Sku): Option[ResourceGroup] =
    sku.getCategory.getResourceGroup match {
      case Cpu.groupName => Some(Cpu)
      case Ram.groupName => Some(Ram)
      case N1Standard.groupName => Some(N1Standard)
      case _ => Option.empty
    }
}
