package cromwell.engine.backend.lsf

import java.nio.file.FileSystem

import cromwell.engine.backend.CallContext
import cromwell.engine.backend.local.LocalCallEngineFunctions

class LsfCallEngineFunctions(fileSystems: List[FileSystem], callContext: CallContext) extends LocalCallEngineFunctions(fileSystems, callContext)
