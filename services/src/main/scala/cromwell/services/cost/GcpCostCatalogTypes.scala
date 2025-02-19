package cromwell.services.cost

import cats.implicits.catsSyntaxValidatedId
import com.google.cloud.billing.v1.Sku
import common.validation.ErrorOr.ErrorOr

import java.util.regex.{Matcher, Pattern}

/*
 * Case class that contains information retrieved from Google about a VM that cromwell has started
 */
case class InstantiatedVmInfo(region: String, machineType: String, preemptible: Boolean)

/*
 * These types reflect hardcoded strings found in a google cost catalog.
 */

sealed trait MachineType { def machineTypeName: String }
case object N1 extends MachineType { override val machineTypeName = "n1" }
case object N2 extends MachineType { override val machineTypeName = "n2" }
case object N2d extends MachineType { override val machineTypeName = "n2d" }

object MachineType {
  def fromSku(sku: Sku): Option[MachineType] = {
    val tokenizedDescription = sku.getDescription.toLowerCase.split(" ")
    if (tokenizedDescription.contains(N1.machineTypeName)) Some(N1)
    else if (tokenizedDescription.contains(N2.machineTypeName)) Some(N2)
    else if (tokenizedDescription.contains(N2d.machineTypeName)) Some(N2d)
    else Option.empty
  }

  // expects a string that looks something like "n1-standard-1" or "custom-1-4096"
  def fromGoogleMachineTypeString(machineTypeString: String): ErrorOr[MachineType] = {
    val mType = machineTypeString.toLowerCase
    if (mType.startsWith("n1-")) N1.validNel
    else if (mType.startsWith("n2d-")) N2d.validNel
    else if (mType.startsWith("n2-")) N2.validNel
    else if (mType.startsWith("custom-")) N1.validNel // by convention
    else s"Unrecognized machine type: $machineTypeString".invalidNel
  }

  def extractCoreCountFromMachineTypeString(machineTypeString: String): ErrorOr[Int] = {
    // Regex to capture second-to-last hyphen-delimited token as number
    val pattern: Pattern = Pattern.compile("-(\\d+)-[^-]+$")
    val matcher: Matcher = pattern.matcher(machineTypeString)
    if (matcher.find()) {
      matcher.group(1).toInt.validNel
    } else {
      s"Could not extract core count from ${machineTypeString}".invalidNel
    }
  }
  def extractRamMbFromMachineTypeString(machineTypeString: String): ErrorOr[Int] = {
    // Regular expression to match the number after a hyphen at the end of the string
    val pattern: Pattern = Pattern.compile("-(\\d+)$")
    val matcher: Matcher = pattern.matcher(machineTypeString);
    if (matcher.find()) {
      matcher.group(1).toInt.validNel
    } else {
      s"Could not extract Ram MB count from ${machineTypeString}".invalidNel
    }
  }
}

sealed trait UsageType { def typeName: String }
case object OnDemand extends UsageType { override val typeName = "ondemand" }
case object Preemptible extends UsageType { override val typeName = "preemptible" }

object UsageType {
  def fromSku(sku: Sku): Option[UsageType] =
    sku.getCategory.getUsageType.toLowerCase match {
      case OnDemand.typeName => Some(OnDemand)
      case Preemptible.typeName => Some(Preemptible)
      case _ => Option.empty
    }
  def fromBoolean(isPreemptible: Boolean): UsageType = isPreemptible match {
    case true => Preemptible
    case false => OnDemand
  }

}

sealed trait MachineCustomization { def customizationName: String }
case object Custom extends MachineCustomization { override val customizationName = "custom" }
case object Predefined extends MachineCustomization { override val customizationName = "predefined" }

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
    val tokenizedDescription = sku.getDescription.toLowerCase.split(" ")

    // ex. "N1 Predefined Instance Core running in Montreal"
    if (tokenizedDescription.contains(Predefined.customizationName)) Some(Predefined)
    // ex. "N2 Custom Instance Core running in Paris"
    else if (tokenizedDescription.contains(Custom.customizationName)) Some(Custom)
    // ex. "N2 Instance Core running in Paris"
    else Some(Predefined)
  }
}

sealed trait ResourceType { def groupName: String }
case object Cpu extends ResourceType { override val groupName = "cpu" }
case object Ram extends ResourceType { override val groupName = "ram" }

object ResourceType {
  def fromSku(sku: Sku): Option[ResourceType] = {
    val tokenizedDescription = sku.getDescription.toLowerCase.split(" ")
    sku.getCategory.getResourceGroup.toLowerCase match {
      case Cpu.groupName => Some(Cpu)
      case Ram.groupName => Some(Ram)
      case "n1standard" if tokenizedDescription.contains("ram") => Some(Ram)
      case "n1standard" if tokenizedDescription.contains("core") => Some(Cpu)
      case _ => Option.empty
    }
  }
}
