package cromwell.backend.impl.bcs

import akka.actor.Actor
import cromwell.backend.standard.StandardCachingActorHelper
import cromwell.core.logging.JobLogging
import cromwell.core.path.Path

trait BcsJobCachingActorHelper extends StandardCachingActorHelper {
  this: Actor with JobLogging =>
  lazy val initializationData: BcsBackendInitializationData = {
    backendInitializationDataAs[BcsBackendInitializationData]
  }

  lazy val bcsClient = initializationData.bcsConfiguration.bcsClient.getOrElse(throw new RuntimeException("no bcs client available"))

  lazy val bcsWorkflowPaths: BcsWorkflowPaths = workflowPaths.asInstanceOf[BcsWorkflowPaths]

  lazy val bcsJobPaths: BcsJobPaths = jobPaths.asInstanceOf[BcsJobPaths]

  lazy val bcsConfiguration: BcsConfiguration = initializationData.bcsConfiguration

  lazy val runtimeAttributes = BcsRuntimeAttributes(validatedRuntimeAttributes, bcsConfiguration.runtimeConfig)

  lazy val callRootPath: Path = bcsJobPaths.callExecutionRoot

  lazy val returnCodeFilename: String = bcsJobPaths.returnCodeFilename
  lazy val returnCodeGcsPath: Path = bcsJobPaths.returnCode
  lazy val standardPaths = bcsJobPaths.standardPaths
  lazy val bcsStdoutFile: Path = standardPaths.output
  lazy val bcsStderrFile: Path = standardPaths.error

  //lazy val bcsCommandLine = "bash -c $(pwd)/cromwell_bcs && sync"
  lazy val bcsCommandLine = "python -u cromwell_bcs.py"
}
