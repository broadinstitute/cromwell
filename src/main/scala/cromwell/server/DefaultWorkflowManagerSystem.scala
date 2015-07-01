package cromwell.server

import cromwell.engine.db.DataAccess

case class DefaultWorkflowManagerSystem() extends WorkflowManagerSystem {
  def dataAccess: DataAccess = DataAccess()
}
