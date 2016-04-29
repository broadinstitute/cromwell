package cromwell.engine.db

@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
case class CallCacheHit(workflowId: String, callName: String)
