package cromwell.backend.sfs

import akka.actor.{ActorSystem, Props}
import akka.testkit.TestActorRef
import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.io.WorkflowPathsWithDocker
import cromwell.backend.validation.{DockerValidation, RuntimeAttributesValidation}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}

import scala.concurrent.Promise

class TestLocalAsyncJobExecutionActor(override val params: SharedFileSystemAsyncJobExecutionActorParams)
  extends BackgroundAsyncJobExecutionActor {
  override lazy val processArgs = {
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
                   (implicit system: ActorSystem): SharedFileSystemJobExecutionActor = {
    createBackendRef(jobDescriptor, configurationDescriptor).underlyingActor
  }

  def createBackendRef(jobDescriptor: BackendJobDescriptor, configurationDescriptor: BackendConfigurationDescriptor)
                      (implicit system: ActorSystem): TestActorRef[SharedFileSystemJobExecutionActor] = {
    val emptyActor = system.actorOf(Props.empty)
    val workflowPaths = new WorkflowPathsWithDocker(jobDescriptor.workflowDescriptor, configurationDescriptor.backendConfig)
    val initializationData = new SharedFileSystemBackendInitializationData(workflowPaths,
      SharedFileSystemValidatedRuntimeAttributesBuilder.default.withValidation(DockerValidation.optional))

    def propsCreator(completionPromise: Promise[BackendJobExecutionResponse]): Props = {
      val params = SharedFileSystemAsyncJobExecutionActorParams(emptyActor, jobDescriptor,
        configurationDescriptor, completionPromise, Option(initializationData))
      Props(classOf[TestLocalAsyncJobExecutionActor], params)
    }

    TestActorRef(new SharedFileSystemJobExecutionActor(
      jobDescriptor, configurationDescriptor, emptyActor, propsCreator))
  }
}
