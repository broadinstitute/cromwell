package cromwell.engine.db

import cromwell.database.obj.Execution

case class ExecutionWithCacheData(execution: Execution, cacheHit: Option[CallCacheHit])
