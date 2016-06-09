package cromwell.engine.workflow.lifecycle

import java.nio.file.Path

import cromwell.backend.BackendWorkflowFinalizationActor.{FinalizationResponse, FinalizationSuccess}
import cromwell.core.PathFactory
import cromwell.engine.{EngineWorkflowDescriptor, WdlFunctions}

import scala.concurrent.{ExecutionContext, Future}

trait EngineWorkflowCopyFinalizationActor extends EngineWorkflowFinalizationActor {

  private lazy val pathFactory = new PathFactory {}
  private lazy val expressionLanguageFunctions = new WdlFunctions(workflowDescriptor.backendDescriptor.workflowOptions)
  private lazy val fileSystems = expressionLanguageFunctions.fileSystems

  protected def workflowDescriptor: EngineWorkflowDescriptor

  protected def copyFiles(): Unit

  final override def afterAll()(implicit ec: ExecutionContext): Future[FinalizationResponse] = Future {
    copyFiles()
    FinalizationSuccess
  }

  protected def getWorkflowOption(key: String): Option[String] = {
    val workflowOptions = workflowDescriptor.backendDescriptor.workflowOptions
    workflowOptions.get(key).toOption
  }

  protected def convertStringToPath(str: String): Path = {
    pathFactory.buildPath(str, fileSystems)
  }
}
