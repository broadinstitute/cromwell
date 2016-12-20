package cromwell.backend.sfs

import akka.actor.{Actor, ActorSystem, Props}
import akka.testkit.TestActorRef
import cromwell.backend.standard._
import better.files.File
import cromwell.backend.io.WorkflowPathsWithDocker
import cromwell.backend.standard.{DefaultStandardSyncExecutionActorParams, StandardAsyncExecutionActorParams, StandardSyncExecutionActor}
import cromwell.backend.validation.{DockerValidation, RuntimeAttributesValidation}
import cromwell.backend.{BackendConfigurationDescriptor, BackendJobDescriptor}
import cromwell.services.io.promise.{ReadAsStringCommandPromise, SizeCommandPromise}

import scala.util.Try

class TestLocalAsyncJobExecutionActor(override val standardParams: StandardAsyncExecutionActorParams)
  extends BackgroundAsyncJobExecutionActor {
  override lazy val processArgs: SharedFileSystemCommand = {
    val script = jobPaths.script.toString
    if (isDockerRun) {
      val docker = RuntimeAttributesValidation.extract(DockerValidation.instance, validatedRuntimeAttributes)
      val cwd = jobPaths.callRoot.toString
      val dockerCwd = jobPathsWithDocker.callDockerRoot.toString
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
    val ioActor = system.actorOf(Props(new TestIoActor))
    val workflowPaths = new WorkflowPathsWithDocker(jobDescriptor.workflowDescriptor, configurationDescriptor.backendConfig)
    val initializationData = new StandardInitializationData(workflowPaths,
      StandardValidatedRuntimeAttributesBuilder.default.withValidation(DockerValidation.optional))
    val asyncClass = classOf[TestLocalAsyncJobExecutionActor]

    val params = DefaultStandardSyncExecutionActorParams(SharedFileSystemAsyncJobExecutionActor.JobIdKey, ioActor, jobDescriptor,
      configurationDescriptor, Option(initializationData), None, asyncClass)

    TestActorRef(new StandardSyncExecutionActor(params))
  }
}

class TestIoActor extends Actor {
  override def receive: Receive = {
    case command: ReadAsStringCommandPromise =>
      command.promise.complete(Try(File(command.file).contentAsString))
      ()
    case command: SizeCommandPromise =>
      command.promise.complete(Try(File(command.file).size))
      ()
  }
}
