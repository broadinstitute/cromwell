package cromwell.backend.impl.local

import java.nio.file.FileSystems

import cromwell.backend.BackendJobExecutionActor.BackendJobExecutionResponse
import cromwell.backend.BackendLifecycleActor.JobAbortResponse
import cromwell.backend._
import cromwell.backend.impl.local.LocalJobExecutionFSM.{Abort, Run}
import org.slf4j.LoggerFactory
import wdl4s._

import scala.concurrent.{Future, Promise}
import scala.language.postfixOps

object LocalJobExecutionActor {
  val ProcessKilledCode = 143
  val logger = LoggerFactory.getLogger("LocalBackend")
  // TODO Support GCS ?
  val fileSystems = List(FileSystems.getDefault)

  private def splitFqn(fullyQualifiedName: FullyQualifiedName): (String, String) = {
    val lastIndex = fullyQualifiedName.lastIndexOf(".")
    (fullyQualifiedName.substring(0, lastIndex), fullyQualifiedName.substring(lastIndex + 1))
  }

  case class Command(argv: Seq[String]) {
    override def toString = argv.map(s => "\"" + s + "\"").mkString(" ")
  }
}

class LocalJobExecutionActor(override val jobDescriptor: BackendJobDescriptor,
                   override val configurationDescriptor: BackendConfigurationDescriptor) extends BackendJobExecutionActor {

  val executionFSM = context.actorOf(LocalJobExecutionFSM.props(jobDescriptor, backendConfiguration))

  override def execute: Future[BackendJobExecutionResponse] = {
    val executePromise = Promise[BackendJobExecutionResponse]()
    executionFSM ! Run(executePromise)
    executePromise.future
  }

  override def abortJob: Future[JobAbortResponse] = {
    val abortPromise = Promise[JobAbortResponse]()
    executionFSM ! Abort(abortPromise)
    abortPromise.future
  }

  override def recover: Future[BackendJobExecutionResponse] = execute
}
