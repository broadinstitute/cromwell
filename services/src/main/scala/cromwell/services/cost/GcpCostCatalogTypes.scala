package cromwell.services.cost

import com.google.cloud.billing.v1.Sku

import java.util.regex.{Matcher, Pattern}
import scala.util.{Failure, Success, Try}

/*
 * Case class that contains information retrieved from Google about a VM that cromwell has started
 */
case class InstantiatedVmInfo(region: String, machineType: String, preemptible: Boolean)
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
    if (machineTypeString.toLowerCase().startsWith("n1-")) Some(N1)
    else if (machineTypeString.toLowerCase().startsWith("n2d-")) Some(N2d)
    else if (machineTypeString.toLowerCase().startsWith("n2-")) Some(N2)
    else if (machineTypeString.toLowerCase().startsWith("custom-")) Some(N1) // by convention
    else {
      println(s"Error: Unrecognized machine type: $machineTypeString")
      None
    }

  def extractCoreCountFromMachineTypeString(machineTypeString: String): Try[Int] = {
    val pattern: Pattern = Pattern.compile("-(\\d+)")
    val matcher: Matcher = pattern.matcher(machineTypeString)
    if (matcher.find()) {
      Success(matcher.group(1).toInt)
    } else {
      Failure(new IllegalArgumentException(s"Could not extract core count from ${machineTypeString}"))
    }
  }
  def extractRamMbFromMachineTypeString(machineTypeString: String): Try[Int] = {
    // Regular expression to match the number after the dash at the end of the string
    // TODO add test
    val pattern: Pattern = Pattern.compile(".*?-(\\d+)$")
    val matcher: Matcher = pattern.matcher(machineTypeString);
    if (matcher.find()) {
      Success(matcher.group(1).toInt)
    } else {
      Failure(new IllegalArgumentException(s"Could not Ram MB count from ${machineTypeString}"))
    }
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
  def fromBoolean(isPreemptible: Boolean): UsageType = isPreemptible match {
    case true => Preemptible
    case false => OnDemand
  }

}
sealed trait UsageType { def typeName: String }
case object OnDemand extends UsageType { override val typeName = "OnDemand" }
case object Preemptible extends UsageType { override val typeName = "Preemptible" }

object MachineCustomization {
  def fromMachineTypeString(machineTypeString: String): MachineCustomization =
    if (machineTypeString.toLowerCase.contains("custom")) Custom else Predefined

  /*
  The cost catalog is annoyingly inconsistent and unstructured in this area.
   - For N1 machines, only predefined SKUs are included, and they have "Predefined" in their description strings.
   We will eventually fall back to using these SKUs for custom machines, but accurately represent them as Predefined here.
   - For non-N1 machines, both custom and predefined SKUs are included, custom ones include "Custom" in their description
   strings and predefined SKUs are only identifiable by the absence of "Custom."
   */
  def fromSku(sku: Sku): Option[MachineCustomization] = {
    val tokenizedDescription = sku.getDescription.split(" ")

    // ex. "N1 Predefined Instance Core running in Montreal"
    if (tokenizedDescription.contains(Predefined.customizationName)) Some(Predefined)
    // ex. "N2 Custom Instance Core running in Paris"
    else if (tokenizedDescription.contains(Custom.customizationName)) Some(Custom)
    // ex. "N2 Instance Core running in Paris"
    else Some(Predefined)
  }
}
sealed trait MachineCustomization { def customizationName: String }
case object Custom extends MachineCustomization { override val customizationName = "Custom" }
case object Predefined extends MachineCustomization { override val customizationName = "Predefined" }

object ResourceType {
  def fromSku(sku: Sku): Option[ResourceType] = {
    val tokenizedDescription = sku.getDescription.split(" ")
    sku.getCategory.getResourceGroup match {
      case Cpu.groupName => Some(Cpu)
      case Ram.groupName => Some(Ram)
      case "N1Standard" if tokenizedDescription.contains("Ram") => Some(Ram)
      case "N1Standard" if tokenizedDescription.contains("Core") => Some(Cpu)
      case _ => Option.empty
    }
  }
}
sealed trait ResourceType { def groupName: String }
case object Cpu extends ResourceType { override val groupName = "CPU" }
case object Ram extends ResourceType { override val groupName = "RAM" }
