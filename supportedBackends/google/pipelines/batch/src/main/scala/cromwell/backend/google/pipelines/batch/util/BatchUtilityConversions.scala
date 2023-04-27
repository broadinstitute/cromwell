package cromwell.backend.google.pipelines.batch.util

import com.google.cloud.batch.v1.AllocationPolicy.{Accelerator, AttachedDisk, Disk, ProvisioningModel}
import com.google.cloud.batch.v1.Volume
import cromwell.backend.google.pipelines.batch.io.{DiskType, GcpBatchAttachedDisk, GcpBatchReferenceFilesDisk}
import cromwell.backend.google.pipelines.batch.models.{GcpBatchConfigurationAttributes, GcpBatchRuntimeAttributes, GpuResource}
import wom.format.MemorySize

trait BatchUtilityConversions {


  // construct zones string
  def toZonesPath(zones: Vector[String]): String = {
    "zones/" + zones.mkString(",")
  }

  // lowercase text to match gcp label requirements
  def toLabel(text: String): String = {
    text.toLowerCase
  }

  // creates batch run time location from zones entered in runtime.  This is needed if network not defined by user to place on right network.
  def toBatchRunLocation(zones: Vector[String]): String = {
    val parts = zones.mkString(",").split("-")
    parts(0) + "-" + parts(1)
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

  def toVolumes(disks: Seq[GcpBatchAttachedDisk]): List[Volume] = disks.map(toVolume).toList

  def toVolume(disk: GcpBatchAttachedDisk): Volume = {
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

  private def toDisk(disk: GcpBatchAttachedDisk): AttachedDisk = {
    val googleDisk = Disk
      .newBuilder
      .setSizeGb(disk.sizeGb.toLong)
      .setType(toBatchDiskType(disk.diskType))

    disk match {
      case refDisk: GcpBatchReferenceFilesDisk =>
        googleDisk.setImage(refDisk.image)
          .build
      case _ =>
        googleDisk.build
    }

    val googleAttachedDisk = AttachedDisk
      .newBuilder
      .setDeviceName(disk.name)
      .setNewDisk(googleDisk)
      .build
    googleAttachedDisk

  }

  private def toBatchDiskType(diskType: DiskType) = diskType match {
    case DiskType.HDD => "pd-standard"
    case DiskType.SSD => "pd-ssd"
    case DiskType.LOCAL => "local-ssd"
  }

  def toVpcNetwork(batchAttributes: GcpBatchConfigurationAttributes): String = {
    batchAttributes.virtualPrivateCloudConfiguration.literalsOption.map { vpcNetworks =>
      vpcNetworks.network
    }.getOrElse(s"projects/${batchAttributes.project}/global/networks/default")
  }

  def toVpcSubnetwork(batchAttributes: GcpBatchConfigurationAttributes, runtimeAttributes: GcpBatchRuntimeAttributes): String = {

    val batchRunLocation = toBatchRunLocation(runtimeAttributes.zones)

    batchAttributes.virtualPrivateCloudConfiguration.labelsOption.map { vpcNetworks =>
      vpcNetworks.subnetwork.getOrElse("default")
    }.getOrElse(s"projects/${batchAttributes.project}/regions/${batchRunLocation}/subnetworks/default")
  }

  def convertGbToMib(runtimeAttributes: GcpBatchRuntimeAttributes): Long = {
    (runtimeAttributes.bootDiskSize * 953.7).toLong
  }

  // Create accelerators for GPUs
  def toAccelerator(gpuResource: GpuResource): Accelerator.Builder = Accelerator.newBuilder.setCount(gpuResource.gpuCount.value.toLong).setType(gpuResource.gpuType.toString)

}
