package cromwell.engine.backend.sge

import cromwell.engine.CallContext
import cromwell.engine.backend.local.LocalCallEngineFunctions
import cromwell.engine.io.IoInterface

class SgeCallEngineFunctions(interface: IoInterface, callContext: CallContext) extends LocalCallEngineFunctions(interface, callContext)
