package cromwell.services.cost

import com.google.cloud.billing.v1.Sku

/*
 * These types reflect hardcoded strings found in a google cost catalog.
 */
object MachineType {
  def fromSku(sku: Sku): Option[MachineType] = {
    val tokenizedDescription = sku.getDescription.split(" ")
    if (tokenizedDescription.contains(N1.machineTypeName)) Some(N1)
    else if (tokenizedDescription.contains(N2.machineTypeName)) Some(N2)
    else if (tokenizedDescription.contains(N2d.machineTypeName)) Some(N2d)
    else Option.empty
  }

  // expects a string that looks something like "n1-standard-1" or "custom-1-4096"
  def fromGoogleMachineTypeString(machineTypeString: String): Option[MachineType] =
    if (machineTypeString.startsWith("n1")) Some(N1)
    else if (machineTypeString.startsWith("n2d")) Some(N2d)
    else if (machineTypeString.startsWith("n2")) Some(N2)
    else if (machineTypeString.startsWith("custom"))
      None // TODO: should this be n1? Make a 'custom' type? Combine with MachineCustomization?
    else {
      println(s"Error: Unrecognized machine type: $machineTypeString")
      None
    }

}
sealed trait MachineType { def machineTypeName: String }
case object N1 extends MachineType { override val machineTypeName = "N1" }
case object N2 extends MachineType { override val machineTypeName = "N2" }
case object N2d extends MachineType { override val machineTypeName = "N2D" }

object UsageType {
  def fromSku(sku: Sku): Option[UsageType] =
    sku.getCategory.getUsageType match {
      case OnDemand.typeName => Some(OnDemand)
      case Preemptible.typeName => Some(Preemptible)
      case _ => Option.empty
    }
}
sealed trait UsageType { def typeName: String }
case object OnDemand extends UsageType { override val typeName = "OnDemand" }
case object Preemptible extends UsageType { override val typeName = "Preemptible" }

object MachineCustomization {
  def fromSku(sku: Sku): Option[MachineCustomization] = {
    val tokenizedDescription = sku.getDescription.split(" ")
    if (tokenizedDescription.contains(Predefined.customizationName)) Some(Predefined)
    else if (tokenizedDescription.contains(Custom.customizationName)) Some(Custom)
    else Option.empty
  }
}
sealed trait MachineCustomization { def customizationName: String }
case object Custom extends MachineCustomization { override val customizationName = "Custom" }
case object Predefined extends MachineCustomization { override val customizationName = "Predefined" }

object ResourceGroup {
  def fromSku(sku: Sku): Option[ResourceGroup] =
    sku.getCategory.getResourceGroup match {
      case Cpu.groupName => Some(Cpu)
      case Ram.groupName => Some(Ram)
      case N1Standard.groupName => Some(N1Standard)
      case _ => Option.empty
    }
}
sealed trait ResourceGroup { def groupName: String }
case object Cpu extends ResourceGroup { override val groupName = "CPU" }
case object Ram extends ResourceGroup { override val groupName = "RAM" }
case object N1Standard extends ResourceGroup { override val groupName = "N1Standard" }
