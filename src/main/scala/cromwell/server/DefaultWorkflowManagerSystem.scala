package cromwell.server

import cromwell.engine.backend.Backend

case class DefaultWorkflowManagerSystem(backend: Backend) extends WorkflowManagerSystem