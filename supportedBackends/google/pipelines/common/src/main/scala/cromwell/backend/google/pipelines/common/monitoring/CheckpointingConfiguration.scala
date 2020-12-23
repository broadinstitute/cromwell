package cromwell.backend.google.pipelines.common.monitoring

import cromwell.backend.BackendJobDescriptor
import cromwell.backend.io.WorkflowPaths
import cromwell.core.path.Path

import scala.concurrent.duration.FiniteDuration

final class CheckpointingConfiguration(jobDescriptor: BackendJobDescriptor,
                                       workflowPaths: WorkflowPaths,
                                       commandDirectory: Path,
                                       checkpointInterval: FiniteDuration
                                      ) {
  def checkpointFileCloud(checkpointFileName: String): String = {
    // The checkpoint file for ANY attempt always goes in the "attempt 1" directory. That way we guarantee that
    // every attempt is able to recover from the single source of checkpointing truth.
    workflowPaths.toJobPaths(jobDescriptor.key.copy(attempt = 1), jobDescriptor.workflowDescriptor)
      .callExecutionRoot.resolve("__checkpointing").resolve(checkpointFileName).toAbsolutePath.pathAsString
  }
  def tmpCheckpointFileCloud(checkpointFileName: String): String = checkpointFileCloud(checkpointFileName) + "-tmp"

  def checkpointFileLocal(checkpointFileName: String): String = {
    commandDirectory.resolve(checkpointFileName).toAbsolutePath.pathAsString
  }
  def tmpCheckpointFileLocal(checkpointFileName: String): String = checkpointFileLocal(checkpointFileName) + "-tmp"

  def localizePreviousCheckpointCommand(checkpointFileName: String): String = {
    val local = checkpointFileLocal(checkpointFileName)
    val cloud = checkpointFileCloud(checkpointFileName)

    s"gsutil cp $cloud $local || touch $local"
  }

  def checkpointingCommand(checkpointFilename: String, multilineActionSquasher: String => String): List[String] = {
    val local = checkpointFileLocal(checkpointFilename)
    val localTmp = tmpCheckpointFileLocal(checkpointFilename)
    val cloud = checkpointFileCloud(checkpointFilename)
    val cloudTmp = tmpCheckpointFileCloud(checkpointFilename)

    val checkpointUploadScript =
      s"""touch $local
         |while true
         |do
         |  # Attempt to make a local copy of the checkpoint file
         |  echo "CHECKPOINTING: Making local copy of $local"
         |  COPY_SUCCESS="false"
         |  while [ "$$COPY_SUCCESS" != "true" ]
         |  do
         |    PRE_COPY_TIMESTAMP="$$(stat -c'%Z' $local)"
         |    cp $local $localTmp
         |    if [ "$$PRE_COPY_TIMESTAMP" == "$$(stat -c'%Z' $local)" ]
         |    then
         |      COPY_SUCCESS="true"
         |    else
         |      echo "CHECKPOINTING: $local was modified while trying to make a local copy. Will retry in 10s..."
         |      sleep 10
         |    fi
         |  done
         |
         |  # Perform the upload:
         |  echo "CHECKPOINTING: Uploading new checkpoint content"
         |  gsutil -m mv $localTmp $cloudTmp
         |  echo "CHECKPOINTING: Replacing cloud checkpoint file with new content"
         |  gsutil -m mv $cloudTmp $cloud
         |  echo "CHECKPOINTING: Sleeping for ${checkpointInterval.toString} before next checkpoint"
         |  sleep ${checkpointInterval.toSeconds}
         |done""".stripMargin

    List(
      "/bin/bash",
      "-c",
      multilineActionSquasher(checkpointUploadScript)
    )
  }
}
