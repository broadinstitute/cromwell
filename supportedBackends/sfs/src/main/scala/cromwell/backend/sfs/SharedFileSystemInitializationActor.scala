package cromwell.backend.sfs

import cromwell.backend.BackendInitializationData
import cromwell.backend.io.WorkflowPaths
import cromwell.backend.standard.{
  StandardExpressionFunctions,
  StandardInitializationActor,
  StandardInitializationActorParams
}
import cromwell.backend.wfs.WorkflowPathBuilder

import scala.concurrent.Future

class SharedFileSystemInitializationActor(standardParams: StandardInitializationActorParams)
    extends StandardInitializationActor(standardParams) {

  override lazy val pathBuilders =
    standardParams.configurationDescriptor.pathBuildersWithDefault(workflowDescriptor.workflowOptions)

  override lazy val workflowPaths: Future[WorkflowPaths] = pathBuilders map {
    WorkflowPathBuilder.workflowPaths(configurationDescriptor, workflowDescriptor, _)
  }

  override lazy val expressionFunctions: Class[_ <: StandardExpressionFunctions] =
    classOf[SharedFileSystemExpressionFunctions]

  override def beforeAll(): Future[Option[BackendInitializationData]] =
    initializationData map { data =>
      publishWorkflowRoot(data.workflowPaths.workflowRoot.pathAsString)
      data.workflowPaths.workflowRoot.createDirectories()
      Option(data)
    }
}
