package cromwell.core.exception
  
import java.io.IOException
  
class LocalizationException(message: String, cause: Option[Throwable]) extends IOException(message) with CromwellGenericException {
  override lazy val exceptionType = "localization"
}
