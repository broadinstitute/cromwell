package cromwell.backend.impl.bcs

import com.aliyuncs.batchcompute.main.v20151111.BatchComputeClient
import com.aliyuncs.batchcompute.pojo.v20151111.TaskDescription
import wom.values._


class BcsJobSpec extends BcsTestUtilSpec {

  behavior of s"BcsJob"

  val mockBcsClient = mock[BatchComputeClient]
  val name = "cromwell"
  val description = name
  val command = "python main.py"
  val packagePath =  mockPathBuilder.build("oss://bcs-test/worker.tar.gz").get
  val mounts = Seq.empty[BcsMount]
  val envs = Map.empty[String, String]

  it should "have correct name and other basic info" in {
    val job = withRuntime()
    job.jobDesc.getName shouldEqual name
    job.jobDesc.getDescription shouldEqual description
    job.jobDesc.getType shouldEqual "DAG"

    val task = taskWithRuntime()
    task.getParameters.getCommand.getCommandLine shouldEqual command
    task.getParameters.getCommand.getPackagePath shouldEqual packagePath.pathAsString

  }

  it should "have correct auto release option" in {
    val runtime = Map("autoReleasejob" -> WomBoolean(false))
    withRuntime(runtime).jobDesc.isAutoRelease shouldBe false
  }

  it should "have correct timeout" in {
    val timeout = 3000
    val runtime = Map("timeout" -> WomInteger(timeout))
    taskWithRuntime(runtime).getTimeout shouldEqual timeout
  }

  it should "have correct mounts" in {
    val src = "oss://bcs-job/dir/"
    val dest = "/home/inputs/"
    val writeSupport = false
    val runtime = Map("mounts" -> WomString(s"$src $dest $writeSupport"))
    taskWithRuntime(runtime).getMounts().getEntries should have size(1)
    taskWithRuntime(runtime).getMounts().getEntries.get(0).getSource shouldBe src
    taskWithRuntime(runtime).getMounts().getEntries.get(0).getDestination shouldBe dest
    taskWithRuntime(runtime).getMounts().getEntries.get(0).isWriteSupport shouldBe writeSupport
  }

  it should "have correct cluster id" in {
    val clusterId = "cls-bcs"
    val runtime = Map("cluster" -> WomString(clusterId))
    taskWithRuntime(runtime).getClusterId shouldEqual clusterId
  }

  it should "have correct docker option" in {
    val dockerImage = "ubuntu/latest"
    val dockerPath = "oss://bcs-reg/ubuntu/"toLowerCase()
    val runtime = Map("dockerTag" -> WomString(s"$dockerImage $dockerPath"))
    taskWithRuntime(runtime).getParameters.getCommand.getEnvVars.get(BcsJob.BcsDockerImageEnvKey) shouldEqual null
    taskWithRuntime(runtime).getParameters.getCommand.getEnvVars.get(BcsJob.BcsDockerPathEnvKey) shouldEqual null
  }

  it should "have correct auto cluster configuration" in {
    val resourceType = "Spot"
    val instanceType = "ecs.c1.large"
    val imageId = "img-centos"
    val spotStrategy = "SpotWithPriceLimit"
    val spotPriceLimit = 0.12
    val cluster = s"$resourceType $instanceType $imageId $spotStrategy $spotPriceLimit"
    val imageIdForCallCaching = "img-ubuntu-vpc"
    val reserveOnFail = true
    val cidr = "172.16.16.0/20"
    val vpcId = "vpc-test"
    val systemDiskType = "cloud"
    val systemDiskSize = 40
    val dataDiskType = "cloud_efficiency"
    val dataDiskSize = 250
    val dataDiskMountPoint = "/home/data/"
    val userDataKey = "key"
    val userDataValue = "value"

    val runtime = Map(
      "cluster" -> WomString(cluster),
      "reserveOnFail" -> WomBoolean(reserveOnFail),
      "vpc" -> WomString(s"$cidr $vpcId"),
      "systemDisk" -> WomString(s"$systemDiskType $systemDiskSize"),
      "dataDisk" -> WomString(s"$dataDiskType $dataDiskSize $dataDiskMountPoint"),
      "userData" -> WomString(s"$userDataKey $userDataValue"),
      "imageId" -> WomString(s"$imageIdForCallCaching")
    )

    val task = taskWithRuntime(runtime)
    a [NullPointerException] should be thrownBy task.getClusterId.isEmpty

    val autoCluster = task.getAutoCluster
    autoCluster.isReserveOnFail shouldEqual reserveOnFail
    autoCluster.getImageId shouldEqual imageIdForCallCaching
    autoCluster.getResourceType shouldEqual resourceType
    autoCluster.getInstanceType shouldEqual instanceType
    autoCluster.getSpotStrategy shouldEqual spotStrategy
    autoCluster.getSpotPriceLimit shouldEqual spotPriceLimit.toFloat

    val vpc = autoCluster.getConfigs.getNetworks.getVpc
    vpc.getVpcId shouldEqual vpcId
    vpc.getCidrBlock shouldEqual cidr

    val systemDisk = autoCluster.getConfigs.getDisks.getSystemDisk
    systemDisk.getType shouldEqual systemDiskType
    systemDisk.getSize shouldEqual systemDiskSize

    val dataDisk = autoCluster.getConfigs.getDisks.getDataDisk
    dataDisk.getType shouldEqual dataDiskType
    dataDisk.getSize shouldEqual dataDiskSize
    dataDisk.getMountPoint shouldEqual dataDiskMountPoint

    val userData = autoCluster.getUserData
    userData.get(userDataKey) shouldEqual userDataValue
  }


  private def withRuntime(runtime: Map[String, WomValue] = Map.empty[String, WomValue]): BcsJob = {
    val runtimeAttributes = createBcsRuntimeAttributes(runtime)
    BcsJob(name, description, command, packagePath, runtimeAttributes.mounts.getOrElse(mounts), envs, runtimeAttributes, None, None, mockBcsClient)
  }

  private def taskWithRuntime(runtime: Map[String, WomValue] = Map.empty[String, WomValue]): TaskDescription = {
    val job = withRuntime(runtime)
    job.jobDesc.getDag.getTasks.get("cromwell")
  }

}
