package cromwell.engine.backend.pbs

import java.nio.file.FileSystem

import cromwell.engine.backend.CallContext
import cromwell.engine.backend.local.LocalCallEngineFunctions

class PbsCallEngineFunctions(fileSystems: List[FileSystem], callContext: CallContext) extends LocalCallEngineFunctions(fileSystems, callContext)
