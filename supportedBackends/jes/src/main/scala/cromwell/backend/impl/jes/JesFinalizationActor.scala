package cromwell.backend.impl.jes

import akka.actor.Props
import better.files._
import cromwell.backend.impl.jes.io._
import cromwell.backend.{BackendWorkflowDescriptor, BackendWorkflowFinalizationActor}
import cromwell.core.retry.Retry
import wdl4s.Call

import scala.concurrent.Future

object JesFinalizationActor {
  def props(workflowDescriptor: BackendWorkflowDescriptor, calls: Seq[Call], jesConfiguration: JesConfiguration) = {
    Props(new JesFinalizationActor(workflowDescriptor, calls, jesConfiguration))
  }
}

class JesFinalizationActor (override val workflowDescriptor: BackendWorkflowDescriptor,
                            override val calls: Seq[Call],
                            jesConfiguration: JesConfiguration) extends BackendWorkflowFinalizationActor {

  override val configurationDescriptor = jesConfiguration.configurationDescriptor

  private val workflowPaths = new JesWorkflowPaths(workflowDescriptor, jesConfiguration)

  override def afterAll(): Future[Unit] = {
    if (jesConfiguration.needAuthFileUpload) deleteAuthenticationFile()
    else Future.successful(())
  }

  private def deleteAuthenticationFile() = {
    val delete = () => Future(workflowPaths.gcsAuthFilePath.delete(false))

    Retry.withRetry(delete, isFatal = isFatalJesException, isTransient = isTransientJesException)(context.system) map { _ => () }
  }
}
