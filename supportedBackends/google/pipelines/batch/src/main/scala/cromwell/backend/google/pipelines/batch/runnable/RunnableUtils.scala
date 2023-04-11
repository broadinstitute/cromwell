package cromwell.backend.google.pipelines.batch.runnable

import com.google.cloud.batch.v1.Runnable

object RunnableUtils {
  /** Start background runnables first, leave the rest as is */
  def sortRunnables(containerSetup: List[Runnable],
                            localization: List[Runnable],
                            userRunnable: List[Runnable],
                            memoryRetryRunnable: List[Runnable],
                            deLocalization: List[Runnable],
                            monitoringSetup: List[Runnable],
                            monitoringShutdown: List[Runnable],
                            checkpointingStart: List[Runnable],
                            checkpointingShutdown: List[Runnable],
                            sshAccess: List[Runnable],
                            isBackground: Runnable => Boolean,
                         ): List[Runnable] = {
    val toBeSortedRunnables = localization ++ userRunnable ++ memoryRetryRunnable ++ deLocalization
    val sortedRunnables = toBeSortedRunnables.sortWith({
      case (runnable, _) => isBackground(runnable)
    })

    sshAccess ++ containerSetup ++ monitoringSetup ++ checkpointingStart ++ sortedRunnables ++ checkpointingShutdown ++ monitoringShutdown
  }
}
