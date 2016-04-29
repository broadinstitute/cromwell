package cromwell.engine.db

import cromwell.database.obj.Execution
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class ExecutionWithCacheData(execution: Execution, cacheHit: Option[CallCacheHit])
