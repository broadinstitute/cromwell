package cromwell.backend.google.pipelines.batch

import com.google.cloud.batch.v1.AllocationPolicy.{Accelerator, ProvisioningModel}

trait BatchUtilityConversions {


  // construct zones string
  def toZonesPath(zones: Vector[String]): String = {
    "zones/" + zones.mkString(",")
  }

  // convert cpu cores to millicores that Batch expects
  def toCpuCores(cpu: Long): Long = {
    cpu * 1000
  }

  // set Standard or Spot instances
  def toProvisioningModel(preemption: Int): ProvisioningModel = preemption compare 0 match {
    case 0 => ProvisioningModel.STANDARD
    case 1 => ProvisioningModel.SPOT
  }


  def toVpcNetwork(batchAttributes: GcpBatchConfigurationAttributes): String = {
    batchAttributes.virtualPrivateCloudConfiguration.labelsOption.map { vpcNetworks =>
        vpcNetworks.network
      }.getOrElse(s"projects/${batchAttributes.project}/global/networks/default")
    }

  def toVpcSubnetwork(batchAttributes: GcpBatchConfigurationAttributes): String = {
    batchAttributes.virtualPrivateCloudConfiguration.labelsOption.map { vpcNetworks =>
      vpcNetworks.subnetwork.getOrElse("default")
    }.getOrElse(s"projects/${batchAttributes.project}/regions/${batchAttributes.location}/subnetworks/default")
  }

  def toBootDiskSizeMb(runtimeAttributes: GcpBatchRuntimeAttributes): Long = {
    (runtimeAttributes.bootDiskSize * 1000).toLong
  }


  // Create accelerators for GPUs
  def toAccelerator(gpuResource: GpuResource): Accelerator.Builder = Accelerator.newBuilder.setCount(gpuResource.gpuCount.value.toLong).setType(gpuResource.gpuType.toString)

}
