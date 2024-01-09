package cromwell.backend.google.batch.runnable

import com.google.cloud.batch.v1.{Runnable, Volume}
import cromwell.backend.google.batch.api.GcpBatchRequestFactory.CreateBatchJobParameters
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.GcsTransferConfiguration

trait MonitoringRunnable {
  def monitoringSetupRunnables(createParameters: CreateBatchJobParameters, volumes: List[Volume])(implicit
    gcsTransferConfiguration: GcsTransferConfiguration
  ): List[Runnable] = {

    val monitoringImageScriptRunnables =
      createParameters.monitoringImage.monitoringImageScriptOption match {
        case Some(script) =>
          val localizeScriptRunnable =
            RunnableBuilder.monitoringImageScriptRunnable(
              script,
              createParameters.monitoringImage.monitoringImageScriptContainerPath,
              volumes
            )
          val describeLocalizeScriptRunnable =
            RunnableBuilder.describeDocker(
              "localizing monitoring image script runnable",
              localizeScriptRunnable
            )
          List(describeLocalizeScriptRunnable, localizeScriptRunnable)
        case None => Nil
      }

    val monitoringImageRunnables =
      createParameters.monitoringImage.monitoringImageOption match {
        case Some(image) =>
          val monitoringImage = image
          val monitoringImageCommand = createParameters.monitoringImage.monitoringImageCommand
          val monitoringImageEnvironment = createParameters.monitoringImage.monitoringImageEnvironment

          val monitoringRunnable = RunnableBuilder.backgroundRunnable(
            monitoringImage,
            monitoringImageCommand,
            monitoringImageEnvironment(volumes.map(_.getMountPath)),
            volumes
          )

          val describeMonitoringRunnable = RunnableBuilder.describeDocker("monitoring runnable", monitoringRunnable)

          List(describeMonitoringRunnable, monitoringRunnable)

        case None => Nil
      }

    (monitoringImageScriptRunnables ++ monitoringImageRunnables).map(_.build)
  }

  def monitoringShutdownRunnables(createParameters: CreateBatchJobParameters): List[Runnable] =
    createParameters.monitoringImage.monitoringImageOption match {
      case Some(_) =>
        val terminationRunnable = RunnableBuilder.terminateBackgroundRunnablesRunnable()

        val describeTerminationRunnable =
          RunnableBuilder.describeDocker("terminate monitoring runnable", terminationRunnable)

        List(describeTerminationRunnable, terminationRunnable).map(_.build)
      case None => Nil
    }
}
