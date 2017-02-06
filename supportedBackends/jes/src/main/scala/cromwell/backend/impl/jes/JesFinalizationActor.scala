package cromwell.backend.impl.jes

import cats.instances.future._
import cats.syntax.functor._
import cromwell.backend._
import cromwell.backend.standard.{StandardFinalizationActor, StandardFinalizationActorParams}
import cromwell.core.CallOutputs
import wdl4s.TaskCall

import scala.concurrent.Future
import scala.language.postfixOps

case class JesFinalizationActorParams
(
  workflowDescriptor: BackendWorkflowDescriptor,
  calls: Set[TaskCall],
  jesConfiguration: JesConfiguration,
  jobExecutionMap: JobExecutionMap,
  workflowOutputs: CallOutputs,
  initializationDataOption: Option[BackendInitializationData]
) extends StandardFinalizationActorParams {
  override val configurationDescriptor: BackendConfigurationDescriptor = jesConfiguration.configurationDescriptor
}

class JesFinalizationActor(val jesParams: JesFinalizationActorParams)
  extends StandardFinalizationActor(jesParams) {

  lazy val jesConfiguration: JesConfiguration = jesParams.jesConfiguration

  override def afterAll(): Future[Unit] = {
    for {
      // NOTE: These are currently in series, not in parallel. Not sure how many threads to throw at finalization
      _ <- deleteAuthenticationFile()
      _ <- super.afterAll()
    } yield ()
  }

  private def deleteAuthenticationFile(): Future[Unit] = {
    (jesConfiguration.needAuthFileUpload, workflowPaths) match {
      case (true, Some(paths: JesWorkflowPaths)) => Future { paths.gcsAuthFilePath.delete() } void
      case _ => Future.successful(())
    }
  }
}
