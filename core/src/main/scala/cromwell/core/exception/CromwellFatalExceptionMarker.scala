package cromwell.core.exception

// Marker trait
trait CromwellFatalExceptionMarker extends CromwellGenericException { this: Throwable => }