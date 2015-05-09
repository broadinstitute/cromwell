package cromwell.engine

import cromwell.binding._

class UnsatisfiedInputsException(diagnostics: Map[FullyQualifiedName, String]) extends RuntimeException
