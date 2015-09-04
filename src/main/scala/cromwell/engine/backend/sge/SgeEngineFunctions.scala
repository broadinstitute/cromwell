package cromwell.engine.backend.sge

import java.nio.file.Path

import cromwell.engine.backend.local.LocalEngineFunctions

class SgeEngineFunctions(cwd: Path, stdout: Path, stderr: Path) extends LocalEngineFunctions(cwd, stdout, stderr)
