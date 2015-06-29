package cromwell.server

import cromwell.engine.db.DataAccess
import cromwell.engine.db.slick.DataAccessController


case class DefaultWorkflowManagerSystem() extends WorkflowManagerSystem {
  def dataAccess: DataAccess = DataAccessController
}