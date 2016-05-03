package cromwell.engine.backend.sge

import java.nio.file.FileSystem

import cromwell.core.OldCallContext
import cromwell.engine.backend.local.OldStyleLocalCallEngineFunctions

class SgeCallEngineFunctionsOld(fileSystems: List[FileSystem], callContext: OldCallContext) extends OldStyleLocalCallEngineFunctions(fileSystems, callContext)
