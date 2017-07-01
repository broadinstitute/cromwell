package cromwell.webservice

import akka.http.scaladsl.server.{PathMatcher0, Route}


object WrappedRoute {

  /**
   * Wraps a route with a prefix.
   * Can optionally serve the wrapped route, and the original unwrapped route, until clients are all switched over.
   *
   * @param unwrappedRoute The original unwrapped route.
   */
  implicit class EnhancedWrappedRoute(val unwrappedRoute: Route) extends AnyVal {
    /**
     *
     * @param wrappedPathPrefix The prefix to wrap unwrapped route. Ex: "api" (implicitly converted to a PathMatcher)
     * @param routeUnwrapped  For legacy reasons, should we also serve up the unwrapped route? Defaults to false.
     * @return The wrappedRoute, followed optionally by the unwrappedRoute if routeUnwrapped is true.
     */
    def wrapped(wrappedPathPrefix: PathMatcher0, routeUnwrapped: Boolean = false): Route = {
      import akka.http.scaladsl.server.directives.PathDirectives._
      import akka.http.scaladsl.server.RouteConcatenation._
      val route = pathPrefix(wrappedPathPrefix) {
        unwrappedRoute
      }
      if (routeUnwrapped) route ~ unwrappedRoute else route
    }
  }

}
