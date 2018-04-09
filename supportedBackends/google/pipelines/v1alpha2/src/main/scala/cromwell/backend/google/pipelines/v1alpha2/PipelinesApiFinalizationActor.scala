package cromwell.backend.google.pipelines.v1alpha2

import cromwell.backend.google.pipelines.common.{PipelinesApiFinalizationActorParams, PipelinesApiWorkflowPaths}
import cromwell.core.io.AsyncIoActorClient

import scala.concurrent.Future

class PipelinesApiFinalizationActor(val jesParams: PipelinesApiFinalizationActorParams)
  extends cromwell.backend.google.pipelines.common.PipelinesApiFinalizationActor(jesParams) with AsyncIoActorClient {

  override def afterAll(): Future[Unit] = {
    for {
      // NOTE: These are currently in series, not in parallel. Not sure how many threads to throw at finalization
      _ <- deleteAuthenticationFile()
      _ <- super.afterAll()
    } yield ()
  }

  private def deleteAuthenticationFile(): Future[Unit] = {
    (jesConfiguration.needAuthFileUpload, workflowPaths) match {
      case (true, Some(paths: PipelinesApiWorkflowPaths)) => asyncIo.deleteAsync(paths.gcsAuthFilePath)
      case _ => Future.successful(())
    }
  }
}
