package cromwell.filesystems.gcs.auth

import com.google.auth.oauth2.GoogleCredentials
import org.scalatest.Assertions._

import scala.util.{Failure, Try}

object GoogleAuthModeSpec {
  def assumeHasApplicationDefaultCredentials(): Unit = {
    tryApplicationDefaultCredentials match {
      case Failure(exception) => cancel(exception.getMessage)
      case _ =>
    }
    ()
  }

  private lazy val tryApplicationDefaultCredentials: Try[Unit] = Try {
    GoogleCredentials.getApplicationDefault
    ()
  }
}
