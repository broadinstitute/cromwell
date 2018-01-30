package cromiam.webservice

import akka.http.scaladsl.marshalling.ToResponseMarshallable
import akka.http.scaladsl.model.HttpRequest
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server._
import cromiam.auth.User

object RequestDirectives {
  /**
    *  FIXME: I see why this was done (as nearly every route uses this) but it becomes awfully limiting to
    *  have it as it means one can't use further directives after this. That probably only matters in a few places,
    *  so this could stay as the "almost always" case. Beyond moving them, leaving these as-is for a future refactor
    */
  def handleRequestWithAuthn(directive: Directive0)(f: (User, HttpRequest) => ToResponseMarshallable): Route = {
    import cromiam.auth.User.requireUser

    directive {
      requireUser { user =>
        toStrictEntity(Timeout) {
          extractRequest { req =>
            complete {
              f.apply(user, req)
            }
          }
        }
      }
    }
  }

  def handlePublicRequest(directive: Directive0)(f: (HttpRequest) => ToResponseMarshallable): Route = {
    directive {
      toStrictEntity(Timeout) {
        extractRequest { req =>
          complete {
            f.apply(req)
          }
        }
      }
    }
  }
}
