package cromwell.backend.impl.bcs

import com.aliyuncs.batchcompute.main.v20151111.BatchComputeClient
import com.aliyuncs.batchcompute.model.v20151111._
import com.aliyuncs.batchcompute.pojo.v20151111._
import cromwell.core.ExecutionEvent
import cromwell.core.path.Path

import collection.JavaConverters._
import scala.util.{Failure, Success, Try}

object BcsJob{
  val BcsDockerImageEnvKey = "BATCH_COMPUTE_DOCKER_IMAGE"
  val BcsDockerPathEnvKey = "BATCH_COMPUTE_DOCKER_REGISTRY_OSS_PATH"
}

final case class BcsJob(name: String,
                  description: String,
                  commandString: String,
                  packagePath: Path,
                  mounts: Seq[BcsMount],
                  envs: Map[String, String],
                  runtime: BcsRuntimeAttributes,
                  stdoutPath: Option[Path],
                  stderrPath: Option[Path],
                  batchCompute: BatchComputeClient) {

  lazy val lazyDisks = new Disks
  lazy val lazyConfigs = new Configs
  lazy val lazyVpc = new VPC
  lazy val lazyTask = new TaskDescription
  lazy val lazyJob = new JobDescription
  lazy val lazyCmd = new Command

  def submit(): Try[String] = Try{
    val request: CreateJobRequest = new CreateJobRequest
    request.setJobDescription(jobDesc)
    val response: CreateJobResponse = batchCompute.createJob(request)
    val jobId = response.getJobId
    jobId
  }

  def getStatus(jobId: String): Try[RunStatus] = Try{
    val request: GetJobRequest = new GetJobRequest
    request.setJobId(jobId)
    val response: GetJobResponse = batchCompute.getJob(request)
    val job = response.getJob
    val status = job.getState
    val message = job.getMessage
    val eventList = Seq[ExecutionEvent]()
    RunStatusFactory.getStatus(jobId, status, Some(message), Some(eventList)) match {
      case Success(status) => status
      case Failure(e) => throw e
    }
  }

  def cancel(jobId: String): Unit = {
    // XXX: Do nothing currently.
  }

  private[bcs] def systemDisk: Option[SystemDisk] = runtime.systemDisk map { disk =>
      val systemDisk = new SystemDisk()
      systemDisk.setType(disk.diskType)
      systemDisk.setSize(disk.sizeInGB)
      systemDisk
  }

  private[bcs] def dataDisk: Option[DataDisk] = runtime.dataDisk map { disk =>
    val dataDisk = new DataDisk
    dataDisk.setType(disk.diskType)
    dataDisk.setSize(disk.sizeInGB)
    dataDisk.setMountPoint(disk.mountPoint)
    dataDisk
  }

  // XXX: maybe more elegant way to reduce two options?
  private[bcs] def disks: Option[Disks] = {
    (systemDisk, dataDisk) match {
      case (Some(sys), Some(data)) =>
        lazyDisks.setSystemDisk(sys)
        lazyDisks.setDataDisk(data)
        Some(lazyDisks)
      case (Some(sys), None) =>
        lazyDisks.setSystemDisk(sys)
        Some(lazyDisks)
      case (None, Some(data)) =>
        lazyDisks.setDataDisk(data)
        Some(lazyDisks)
      case (None, None) => None
    }
  }

  private[bcs] def vpc: Option[VPC] = {
    (runtime.vpc flatMap {v => v.cidrBlock}, runtime.vpc flatMap {v => v.vpcId}) match {
      case (Some(cidr), Some(id)) =>
        lazyVpc.setCidrBlock(cidr)
        lazyVpc.setVpcId(id)
        Some(lazyVpc)
      case (Some(cidr), None) =>
        lazyVpc.setCidrBlock(cidr)
        Some(lazyVpc)
      case (None, Some(vpc)) =>
        lazyVpc.setVpcId(vpc)
        Some(lazyVpc)
      case (None, None) => None
    }
  }

  private[bcs] def configs: Option[Configs] = {
    (vpc, disks) match {
      case (Some(bcsVpc), Some(bcsDisks)) =>
        lazyConfigs.setDisks(bcsDisks)
        val networks = new Networks
        networks.setVpc(bcsVpc)
        lazyConfigs.setNetworks(networks)
        Some(lazyConfigs)
      case (Some(bcsVpc), None) =>
        val networks = new Networks
        networks.setVpc(bcsVpc)
        lazyConfigs.setNetworks(networks)
        Some(lazyConfigs)
      case (None, Some(bcsDisks)) =>
        lazyConfigs.setDisks(bcsDisks)
        Some(lazyConfigs)
      case (None, None) => None
    }
  }

  private[bcs] def params: Parameters = {
    val parames = new Parameters
    lazyCmd.setPackagePath(packagePath.pathAsString)
    lazyCmd.setEnvVars(environments.asJava)
    lazyCmd.setCommandLine(commandString)

    dockers foreach {docker => lazyCmd.setDocker(docker)}
    stdoutPath foreach {path => parames.setStdoutRedirectPath(path.normalize().pathAsString + "/")}
    stderrPath foreach {path => parames.setStderrRedirectPath(path.normalize().pathAsString + "/")}

    parames.setCommand(lazyCmd)
    parames
  }

  private[bcs] def environments: Map[String, String] = {
    runtime.docker match {
      case None =>
        runtime.dockerTag match {
          case Some(docker: BcsDockerWithoutPath) => envs + (BcsJob.BcsDockerImageEnvKey -> docker.image)
          case Some(docker: BcsDockerWithPath) => envs + (BcsJob.BcsDockerPathEnvKey -> docker.path) + (BcsJob.BcsDockerImageEnvKey -> docker.image)
          case _ => envs
        }
      case _ => envs
    }
  }

  val dockers: Option[Command.Docker] = {
    runtime.docker match {
      case Some(docker: BcsDockerWithoutPath) =>
        val dockers = new Command.Docker
        dockers.setImage(docker.image)
        Some(dockers)
      case _ => None
    }
  }

  private[bcs] def jobDesc: JobDescription = {
    lazyJob.setName(name)
    lazyJob.setDescription(description)
    lazyJob.setType("DAG")

    val dag = new DAG
    dag.addTask("cromwell", taskDesc)
    lazyJob.setDag(dag)

    // NOTE: Do NOT set auto release here or we will not be able to get status after the job completes.
    lazyJob.setAutoRelease(false)

    lazyJob
  }

  private[bcs] def taskDesc: TaskDescription = {
    lazyTask.setParameters(params)
    lazyTask.setInstanceCount(1)

    runtime.timeout foreach {timeout => lazyTask.setTimeout(timeout.toLong)}

    val cluster = runtime.cluster getOrElse(throw new IllegalArgumentException("cluster id or auto cluster configuration is mandatory"))
    cluster.fold(handleClusterId, handleAutoCluster)

    val mnts = new Mounts
    mounts foreach  {
      case input: BcsInputMount =>
        mnts.addEntries(input.toBcsMountEntry)
      case output: BcsOutputMount =>
        var srcStr = BcsMount.toString(output.src)
        if (BcsMount.toString(output.dest).endsWith("/") && !srcStr.endsWith("/")) {
          srcStr += "/"
        }
        lazyTask.addOutputMapping(srcStr, BcsMount.toString(output.dest))
    }

    lazyTask.setMounts(mnts)

    lazyTask
  }

  private def handleAutoCluster(config: AutoClusterConfiguration): Unit = {
    val autoCluster = new AutoCluster
    autoCluster.setImageId(runtime.imageId.getOrElse(config.imageId))
    autoCluster.setInstanceType(config.instanceType)
    autoCluster.setResourceType(config.resourceType)

    config.spotStrategy foreach {strategy => autoCluster.setSpotStrategy(strategy)}
    config.spotPriceLimit foreach {priceLimit => autoCluster.setSpotPriceLimit(priceLimit)}
    config.clusterId foreach {clusterId => autoCluster.setClusterId(clusterId)}
    runtime.reserveOnFail foreach {reserve => autoCluster.setReserveOnFail(reserve)}
    val userData = runtime.userData map {datas => Map(datas map {data => data.key -> data.value}: _*)}
    userData foreach {datas => autoCluster.setUserData(datas.asJava)}

    configs foreach (bcsConfigs => autoCluster.setConfigs(bcsConfigs))
    lazyTask.setAutoCluster(autoCluster)
  }

  private def handleClusterId(clusterId: String): Unit = lazyTask.setClusterId(clusterId)
}
