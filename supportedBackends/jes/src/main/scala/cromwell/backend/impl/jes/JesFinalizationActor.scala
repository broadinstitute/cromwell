package cromwell.backend.impl.jes

import akka.actor.ActorRef
import cromwell.backend._
import cromwell.backend.standard.{StandardFinalizationActor, StandardFinalizationActorParams}
import cromwell.core.CallOutputs
import cromwell.core.io.AsyncIo
import cromwell.filesystems.gcs.batch.GcsBatchCommandBuilder
import wdl4s.wdl.WdlTaskCall

import scala.concurrent.Future

case class JesFinalizationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  ioActor: ActorRef,
  calls: Set[WdlTaskCall],
  jesConfiguration: JesConfiguration,
  jobExecutionMap: JobExecutionMap,
  workflowOutputs: CallOutputs,
  initializationDataOption: Option[BackendInitializationData]
) extends StandardFinalizationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

class JesFinalizationActor(val jesParams: JesFinalizationActorParams)
  extends StandardFinalizationActor(jesParams) with AsyncIo with GcsBatchCommandBuilder {

  lazy val jesConfiguration: JesConfiguration = jesParams.jesConfiguration
  
  override def receive = ioReceive orElse super.receive

  override def afterAll(): Future[Unit] = {
    for {
      // NOTE: These are currently in series, not in parallel. Not sure how many threads to throw at finalization
      _ <- deleteAuthenticationFile()
      _ <- super.afterAll()
    } yield ()
  }

  private def deleteAuthenticationFile(): Future[Unit] = {
    (jesConfiguration.needAuthFileUpload, workflowPaths) match {
      case (true, Some(paths: JesWorkflowPaths)) => deleteAsync(paths.gcsAuthFilePath)
      case _ => Future.successful(())
    }
  }

  override def ioActor: ActorRef = jesParams.ioActor
}
