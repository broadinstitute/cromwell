package cromwell.engine.backend.sge

import java.nio.file.FileSystem

import cromwell.engine.backend.CallContext
import cromwell.engine.backend.local.LocalCallEngineFunctions

class SgeCallEngineFunctions(fileSystems: List[FileSystem], callContext: CallContext) extends LocalCallEngineFunctions(fileSystems, callContext)
