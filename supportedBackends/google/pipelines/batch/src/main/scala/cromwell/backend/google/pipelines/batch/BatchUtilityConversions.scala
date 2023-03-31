package cromwell.backend.google.pipelines.batch

import com.google.cloud.batch.v1.AllocationPolicy.{Accelerator, ProvisioningModel, Disk, AttachedDisk}
import com.google.cloud.batch.v1.Volume
import cromwell.backend.google.pipelines.batch.io.{DiskType, GcpBatchAttachedDisk, GcpBatchReferenceFilesDisk}
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

  def toDisks(disks: Seq[GcpBatchAttachedDisk]): List[AttachedDisk] = disks.map(toDisk).toList

  def toMount(disk: GcpBatchAttachedDisk): Volume = {
    val volume = Volume
      .newBuilder
      .setDeviceName(disk.name)
      .setMountPath(disk.mountPoint.pathAsString)


    disk match {
      case _: GcpBatchReferenceFilesDisk =>
        volume
        .addMountOptions("async, rw")
          .build
      case _ =>
        volume
          .build
    }
  }

  def toDisk(disk: GcpBatchAttachedDisk): AttachedDisk = {
    val googleDisk = Disk
      .newBuilder
      .setSizeGb(disk.sizeGb.toLong)
      .setType(toBatchDiskType(disk.diskType))
      .build

    val googleAttachedDisk = AttachedDisk
      .newBuilder
      .setDeviceName(disk.name)
      .setNewDisk(googleDisk)
      .build
    googleAttachedDisk
    //disk match {
    //  case refDisk: GcpBatchReferenceFilesDisk =>
    //    googleDisk.setSourceImage(refDisk.image)
    //  case _ =>
    //    googleDisk
    //}
  }

  private def toBatchDiskType(diskType: DiskType) = diskType match {
    case DiskType.HDD => "pd-standard"
    case DiskType.SSD => "pd-ssd"
    case DiskType.LOCAL => "local-ssd"
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

  def convertGbToMib(runtimeAttributes: GcpBatchRuntimeAttributes): Long = {
    (runtimeAttributes.bootDiskSize * 953.7).toLong
  }


  // Create accelerators for GPUs
  def toAccelerator(gpuResource: GpuResource): Accelerator.Builder = Accelerator.newBuilder.setCount(gpuResource.gpuCount.value.toLong).setType(gpuResource.gpuType.toString)

}
