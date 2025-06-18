package cromwell.backend.google.batch.runnable

import com.google.cloud.batch.v1.{Runnable, Volume}
import common.util.StringUtil._
import cromwell.backend.google.batch.api.GcpBatchRequestFactory.CreateBatchJobParameters
import cromwell.backend.google.batch.models.GcpBatchConfigurationAttributes.GcsTransferConfiguration
import cromwell.backend.google.batch.models.GcpBatchJobPaths.GcsDelocalizationScriptName
import cromwell.backend.google.batch.util.GcpBatchParameterConversions._
import cromwell.backend.google.batch.util.RuntimeOutputMapping
import cromwell.backend.google.batch.util.ToParameter.ops._
import cromwell.core.path.Path
import wom.runtime.WomOutputRuntimeExtractor

import java.util.UUID

trait Delocalization {

  import RunnableBuilder._
  import RunnableCommands._
  import RunnableLabels._
  import RunnableUtils._

  private def runtimeOutputExtractorRunnable(containerCallRoot: String,
                                             outputFile: String,
                                             womOutputRuntimeExtractor: WomOutputRuntimeExtractor
  ): Runnable.Builder = {
    val commands = List(
      "-c",
      // Create the directory where the fofn will be written
      s"mkdir -p $$(dirname $outputFile) && " +
        s"cd $containerCallRoot && " +
        """echo "Runtime output files to be delocalized:" && """ +
        s"${womOutputRuntimeExtractor.command} | tee $outputFile"
    )

    RunnableBuilder
      .withImage(womOutputRuntimeExtractor.dockerImage.getOrElse(CloudSdkImage))
      .withCommand(commands: _*)
      .withEntrypointCommand("/bin/bash")
      .withLabels(Map(Key.Tag -> Value.Delocalization))
  }

  private def delocalizeRuntimeOutputsScript(fofnPath: String, workflowRoot: Path, cloudCallRoot: Path)(implicit
    gcsTransferConfiguration: GcsTransferConfiguration
  ) = {
    val gsutilCommand: String => String = { flag =>
      s"""rm -f $$HOME/.config/gcloud/gce && gsutil -m $flag cp -r $$line "${cloudCallRoot.pathAsString.ensureSlashed}$$gcs_path""""
    }

    def sedStripPrefix(prefix: String) = s"""sed -e "s/^${prefix.ensureSedEscaped}//""""

    // See RuntimeOutputMapping.prefixFilters for more info on why this is needed
    val prefixFilters = RuntimeOutputMapping
      .prefixFilters(workflowRoot)
      .map(sedStripPrefix)
      .mkString(" | ")

    /*
     * Delocalize all the files returned by the runtime output extractor
     */
    s"""|#!/usr/bin/env bash
        |
        |set -x
        |
        |if [ -f $fofnPath ]; then
        |  while IFS= read line
        |  do
        |    gcs_path=$$(echo $$line | $prefixFilters)
        |    (
        |       ${retry(recoverRequesterPaysError(cloudCallRoot)(gsutilCommand))}
        |    )
        |  done  <$fofnPath
        |fi""".stripMargin
  }

  private def delocalizeRuntimeOutputsRunnable(cloudCallRoot: Path,
                                               inputFile: String,
                                               workflowRoot: Path,
                                               volumes: List[Volume]
  )(implicit gcsTransferConfiguration: GcsTransferConfiguration): Runnable.Builder = {
    val command = multiLineCommand(delocalizeRuntimeOutputsScript(inputFile, workflowRoot, cloudCallRoot))
    RunnableBuilder.cloudSdkShellRunnable(command)(volumes = volumes,
                                                   labels = Map(Key.Tag -> Value.Delocalization),
                                                   flags = List.empty
    )
  }

  def deLocalizeRunnables(createParameters: CreateBatchJobParameters, volumes: List[Volume])(implicit
    gcsTransferConfiguration: GcsTransferConfiguration
  ): List[Runnable] = {
    val cloudCallRoot = createParameters.cloudCallRoot
    val callExecutionContainerRoot = createParameters.commandScriptContainerPath.parent

    /*
     * Ideally temporaryFofnForRuntimeOutputFiles should be somewhere else than the execution directory (we could mount anther directory)
     * However because it runs after everything else there's no risk of polluting the task's results and the random ID ensures we don't override anything
     */
    val temporaryFofnDirectoryForRuntimeOutputFiles =
      callExecutionContainerRoot.pathAsString.ensureSlashed + UUID.randomUUID().toString.split("-")(0)
    val temporaryFofnForRuntimeOutputFiles = temporaryFofnDirectoryForRuntimeOutputFiles + "/runtime_output_files.txt"

    val runtimeExtractionRunnables = createParameters.womOutputRuntimeExtractor.toList flatMap { extractor =>
      List(
        runtimeOutputExtractorRunnable(callExecutionContainerRoot.pathAsString,
                                       temporaryFofnForRuntimeOutputFiles,
                                       extractor
        ),
        delocalizeRuntimeOutputsRunnable(cloudCallRoot,
                                         temporaryFofnForRuntimeOutputFiles,
                                         createParameters.cloudWorkflowRoot,
                                         volumes
        )
      )
    }

    val gcsDelocalizationContainerPath =
      createParameters.commandScriptContainerPath.sibling(GcsDelocalizationScriptName)

    val delocalizationLabel = Map(Key.Tag -> Value.Delocalization)
    val runGcsDelocalizationScript = cloudSdkShellRunnable(s"/bin/bash $gcsDelocalizationContainerPath")(
      volumes = volumes,
      labels = delocalizationLabel,
      flags = List.empty
    )

    val annotatedRunnables: List[Runnable.Builder] = runGcsDelocalizationScript ::
      createParameters.outputParameters.flatMap(_.toRunnables(volumes)) ++
      runtimeExtractionRunnables

    // NOTE: papiv2 delocalizes logs from /google but such logs are not available on batch
    // See: https://cloud.google.com/life-sciences/docs/reference/rpc/google.cloud.lifesciences.v2beta
    val all = RunnableBuilder.annotateTimestampedRunnable("delocalization", Value.Delocalization)(
      volumes,
      annotatedRunnables
    )

    all.map(_.build)
  }
}
