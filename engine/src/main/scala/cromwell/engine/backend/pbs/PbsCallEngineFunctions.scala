package cromwell.engine.backend.pbs

import cromwell.engine.CallContext
import cromwell.engine.backend.local.LocalCallEngineFunctions
import cromwell.engine.io.IoInterface

class PbsCallEngineFunctions(interface: IoInterface, callContext: CallContext) extends LocalCallEngineFunctions(interface, callContext)
