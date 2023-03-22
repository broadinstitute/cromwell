package cromwell.backend.google.pipelines.batch

import com.google.cloud.batch.v1.AllocationPolicy.{Accelerator, ProvisioningModel}
import wom.format.MemorySize

trait BatchUtilityConversions {


  // construct zones string
  def toZonesPath(zones: Vector[String]): String = {
    "zones/" + zones.mkString(",")
  }

  // convert cpu cores to millicores that Batch expects
  def toCpuCores(cpu: Long): Long = {
    cpu * 1000
  }

  // convert memory to MiB that Batch expects
  def toMemMib(memory: MemorySize): Long = {
    (memory.amount * 1024).toLong
  }

  // set Standard or Spot instances
  def toProvisioningModel(preemption: Int): ProvisioningModel = preemption compare 0 match {
    case 0 => ProvisioningModel.STANDARD
    case 1 => ProvisioningModel.SPOT
  }

  // Create accelerators for GPUs
  def toAccelerator(gpuResource: GpuResource): Accelerator.Builder = Accelerator.newBuilder.setCount(gpuResource.gpuCount.value.toLong).setType(gpuResource.gpuType.toString)

}
