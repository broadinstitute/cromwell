package cromwell.backend.google.pipelines.common

import cromwell.backend.google.pipelines.common.api.RunStatus
import cromwell.backend.standard.costestimation.CostPollingHelper
import cromwell.services.cost.{Cpu, Custom, MachineType, OnDemand}
import cromwell.services.metadata.CallMetadataKeys

import java.time.OffsetDateTime

class PapiCostPollingHelper(tellMetadataFn: Map[String, Any] => Unit) extends CostPollingHelper[RunStatus] {

  override def extractStartTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList.collectFirst {
      case event if event.name == CallMetadataKeys.VmStartTime => event.offsetDateTime
    }

  override def extractEndTimeFromRunState(pollStatus: RunStatus): Option[OffsetDateTime] =
    pollStatus.eventList.collectFirst {
      case event if event.name == CallMetadataKeys.VmEndTime => event.offsetDateTime
    }

  override def extractVmCostPerHourFromRunState(pollStatus: RunStatus): Option[BigDecimal] =
    pollStatus.instantiatedVmInfo.map { vmInfo =>
      val machineType = MachineType.fromGoogleMachineTypeString(vmInfo.machineType)
      val usageType = OnDemand //TODO: Account for preemptible here
      val machineCustomization = Custom //TODO, also account for predefined
      val resourceGroup = Cpu //TODO, also do RAM. For some reason N1Standard is also a resource group. How to account for that?
      val region = vmInfo.region
      // TODO: Use cost catalog service here. It should take ^ and calculate CPU + RAM cost/hr
      // Failure here is fatal. We have learned all we can from Google, so subsequent attempts wont fair any better.
      3.50
    }

  override def tellMetadata(metadata: Map[String, Any]): Unit = tellMetadataFn(metadata)

}
