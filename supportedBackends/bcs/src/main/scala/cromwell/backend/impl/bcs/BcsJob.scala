package cromwell.backend.impl.bcs

import com.aliyuncs.batchcompute.main.v20151111.BatchComputeClient
import com.aliyuncs.batchcompute.model.v20151111._
import com.aliyuncs.batchcompute.pojo.v20151111._
import cromwell.core.ExecutionEvent
import cromwell.core.path.Path

import collection.JavaConverters._
import scala.util.Try

object BcsJob{
  val BcsDockerImageEnvKey = "BATCH_COMPUTE_DOCKER_IMAGE"
  val BcsDockerPathEnvKey = "BATCH_COMPUTE_DOCKER_REGISTRY_OSS_PATH"
}

case class BcsJob(name: String,
                  description: String,
                  commandString: String,
                  packagePath: Path,
                  mounts: Seq[BcsMount],
                  envs: Map[String, String],
                  runtime: BcsRuntimeAttributes,
                  batchCompute: BatchComputeClient) {
  val jobDesc = new JobDescription
  val task = new TaskDescription

  jobDesc.setName(name)
  jobDesc.setDescription(description)
  jobDesc.setType("DAG")

  val cmd = new Command
  cmd.setCommandLine(commandString)

  var updatedEnvs = runtime.dockerImage match {
    case Some(image) => envs + (BcsJob.BcsDockerImageEnvKey -> image)
    case None => envs
  }

  updatedEnvs = runtime.dockerPath match {
    case Some(path) => updatedEnvs + (BcsJob.BcsDockerPathEnvKey -> path)
    case None => updatedEnvs
  }

  cmd.setPackagePath(packagePath.pathAsString)
  cmd.setEnvVars(updatedEnvs.asJava)

  val params = new Parameters
  params.setCommand(cmd)
  task.setParameters(params)
  task.setInstanceCount(1)

  runtime.clusterId match {
    case Some(id) => task.setClusterId(id)
    case None =>
      val autoOption =  for {
        resourceType <- runtime.resourceType
        instanceType <- runtime.instanceType
        imageId <- runtime.imageId
      } yield {
        val autoCluster = new AutoCluster
        autoCluster.setImageId(imageId)
        autoCluster.setResourceType(resourceType)
        autoCluster.setInstanceType(instanceType)

        runtime.reserveOnFail match {
          case Some(reserve) => autoCluster.setReserveOnFail(reserve)
          case None => autoCluster.setReserveOnFail(BcsRuntimeAttributes.ReserveOnFailDefault)
        }

        var userData: Map[String, String] = Map()
        runtime.userData match {
          case Some(datas) => {
            datas map {
              data =>
                {
                  userData += { data.key -> data.value}
                }
            }
          }
          case None => println("no user data")
        }

        autoCluster.setUserData(userData.asJava)
        autoCluster
      }
      task.setAutoCluster(autoOption.getOrElse(throw new RuntimeException("invalid auto cluster parameters")))
  }

  runtime.timeout match {
    case Some(timeout) => task.setTimeout(timeout.toLong)
    case None => task.setTimeout(21600)
  }

  val dag = new DAG
  dag.addTask("cromwell", task)
  jobDesc.setDag(dag)

  mounts map { mount =>
    mount match {
      case input: BcsInputMount =>
        var destStr = input.dest.pathAsString
        if (input.src.pathAsString.endsWith("/") && !destStr.endsWith("/")) {
          destStr += "/"
        }
        task.addInputMapping(input.src.pathAsString, destStr)
      case output: BcsOutputMount =>
        var srcStr = output.src.pathAsString
        if (output.dest.pathAsString.endsWith("/") && !srcStr.endsWith("/")) {
          srcStr += "/"
        }
        task.addOutputMapping(srcStr, output.dest.pathAsString)
    }
  }

  def submit(): Try[String] = Try{
    val request: CreateJobRequest = new CreateJobRequest
    request.setJobDescription(jobDesc)
    val response: CreateJobResponse = batchCompute.createJob(request)
    val jobId = response.getJobId
    jobId
  }

  def getStatus(jobId: String): Try[RunStatus] = {
    val request: GetJobRequest = new GetJobRequest
    request.setJobId(jobId)
    val response: GetJobResponse = batchCompute.getJob(request)
    val job = response.getJob
    val status = job.getState
    val message = job.getMessage
    val eventList = Seq[ExecutionEvent]()
    RunStatusFactory.getStatus(jobId, status, Some(message), Some(eventList))
  }

  def cancel(jobId: String): Unit = {
    // XXX: Do nothing currently.
  }
}
