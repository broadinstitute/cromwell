package cromiam.auth

import akka.http.scaladsl.model.headers.Authorization
import akka.http.scaladsl.server.Directive1
import org.broadinstitute.dsde.workbench.model.WorkbenchUserId

/**
  * Wraps the concept of an authenticated workbench user including their numeric ID as well as their bearer token
  */
final case class User(userId: WorkbenchUserId, authorization: Authorization)

object User {
  import akka.http.scaladsl.server.Directives._

  /**
    * Obtain both the user id header from the proxy as well as the bearer token and pass that back
    * into the route logic as a User object
    */
  def requireUser: Directive1[User] = {
    (headerValueByName("OIDC_CLAIM_user_id") & headerValuePF { case a: Authorization => a }) tmap { case (userId, auth) =>
      User(WorkbenchUserId(userId), auth)
    }
  }
}
