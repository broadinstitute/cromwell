package cromwell.engine.backend.sge

import java.nio.file.Path

import cromwell.binding.IOInterface
import cromwell.engine.backend.local.LocalEngineFunctions

class SgeEngineFunctions(cwd: Path, stdout: Path, stderr: Path, interface: IOInterface) extends LocalEngineFunctions(cwd, stdout, stderr, interface)
