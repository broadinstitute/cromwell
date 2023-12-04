package cromiam.auth

import akka.http.scaladsl.model.headers.Authorization
import org.broadinstitute.dsde.workbench.model.WorkbenchUserId

/**
  * Wraps the concept of an authenticated workbench user including their numeric ID as well as their bearer token
  */
final case class User(userId: WorkbenchUserId, authorization: Authorization)
