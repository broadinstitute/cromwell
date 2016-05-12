package cromwell.engine.backend.sge

import java.nio.file.FileSystem

import cromwell.core.{OldCallContext, CallContext}
import cromwell.engine.backend.local.OldStyleLocalCallEngineFunctions
@deprecated(message = "This class will not be part of the PBE universe", since = "May 2nd 2016")
class SgeCallEngineFunctions(fileSystems: List[FileSystem], callContext: OldCallContext) extends OldStyleLocalCallEngineFunctions(fileSystems, callContext)
