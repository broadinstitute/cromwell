package cromwell.engine.backend

import cromwell.binding.values.WdlFile
import org.joda.time.DateTime

case class CallMetadata(inputs: Map[String, String],
                        start: DateTime,
                        end: DateTime,
                        jobid: String,
                        rc: Int,
                        stdout: WdlFile,
                        stderr: WdlFile)
