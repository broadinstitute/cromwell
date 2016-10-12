package cromwell.core.path

import cromwell.core.CromwellFatalException

case class PathParsingException(message: String) extends CromwellFatalException(new IllegalArgumentException(message))
