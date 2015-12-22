package cromwell.engine.backend.sge

import java.nio.file.Path

import cromwell.binding.IoInterface
import cromwell.engine.CallContext
import cromwell.engine.backend.local.LocalCallEngineFunctions

class SgeEngineFunctions(cwd: Path, stdout: Path, stderr: Path, interface: IoInterface) extends LocalCallEngineFunctions(interface, new CallContext(cwd.toString, stdout.toString, stderr.toString))
