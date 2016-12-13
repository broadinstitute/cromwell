package cromwell.backend.sfs

import akka.actor.{ActorSystem, Props}
import akka.testkit.TestActorRef
import cromwell.backend.standard.{StandardAsyncExecutionActorParams, StandardSyncExecutionActor, DefaultStandardSyncExecutionActorParams}
import cromwell.backend.io.WorkflowPathsWithDocker
import cromwell.backend.validation.{DockerValidation, RuntimeAttributesValidation}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}

class TestLocalAsyncJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackgroundAsyncJobExecutionActor {
  override lazy val processArgs: SharedFileSystemCommand = {
    val script = jobPaths.script.toString
    if (isDockerRun) {
      val docker = RuntimeAttributesValidation.extract(DockerValidation.instance, validatedRuntimeAttributes)
      val cwd = jobPaths.callRoot.toString
      val dockerCwd = jobPaths.callDockerRoot.toString
      SharedFileSystemCommand("/bin/bash", "-c",
        s"docker run --rm -v $cwd:$dockerCwd -i $docker /bin/bash < $script")
    } else {
      SharedFileSystemCommand("/bin/bash", script)
    }
  }
}

object TestLocalAsyncJobExecutionActor {
  def createBackend(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor)
                   (implicit system: ActorSystem): StandardSyncExecutionActor = {
    createBackendRef(jobDescriptor, configurationDescriptor).underlyingActor
  }

  def createBackendRef(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor)
                      (implicit system: ActorSystem): TestActorRef[StandardSyncExecutionActor] = {
    val emptyActor = system.actorOf(Props.empty)
    val workflowPaths = new WorkflowPathsWithDocker(jobDescriptor.workflowDescriptor, configurationDescriptor.backendConfig)
    val initializationData = new SharedFileSystemBackendInitializationData(workflowPaths,
      SharedFileSystemValidatedRuntimeAttributesBuilder.default.withValidation(DockerValidation.optional))
    val asyncClass = classOf[TestLocalAsyncJobExecutionActor]

    val params = DefaultStandardSyncExecutionActorParams(SharedFileSystemJob.JobIdKey, emptyActor, jobDescriptor,
      configurationDescriptor, Option(initializationData), asyncClass)

    TestActorRef(new StandardSyncExecutionActor(params))
  }
}
