package cromwell.webservice

import akka.http.scaladsl.model.StatusCodes
import akka.http.scaladsl.server.Route
import com.typesafe.config.Config
import net.ceedubs.ficus.Ficus._
import akka.http.scaladsl.server.Directives._

/**
 * Serves up the swagger UI from org.webjars/swagger-ui.
 */
trait SwaggerUiHttpService {
  /**
   * @return The version of the org.webjars/swagger-ui artifact. For example "2.1.1".
   */
  def swaggerUiVersion: String

  /**
   * Informs the swagger UI of the base of the application url, as hosted on the server.
   * If your entire app is served under "http://myserver/myapp", then the base URL is "/myapp".
   * If the app is served at the root of the application, leave this value as the empty string.
   *
   * @return The base URL used by the application, or the empty string if there is no base URL. For example "/myapp".
   */
  def swaggerUiBaseUrl: String = ""

  /**
   * @return The path to the swagger UI html documents. For example "swagger"
   */
  def swaggerUiPath: String = "swagger"

  /**
   * The path to the actual swagger documentation in either yaml or json, to be rendered by the swagger UI html.
   *
   * @return The path to the api documentation to render in the swagger UI.
   *         For example "api-docs" or "swagger/lenthall.yaml".
   */
  def swaggerUiDocsPath: String = "api-docs"

  /**
   * @return When true, if someone requests / (or /baseUrl if setup), redirect to the swagger UI.
   */
  def swaggerUiFromRoot: Boolean = true

  private def routeFromRoot: Route = get {
    pathEndOrSingleSlash {
      // Redirect / to the swagger UI
      redirect(s"$swaggerUiBaseUrl/$swaggerUiPath", StatusCodes.TemporaryRedirect)
    }
  }

  /**
   * Serves up the swagger UI only. Redirects requests to the root of the UI path to the index.html.
   *
   * @return Route serving the swagger UI.
   */
  final def swaggerUiRoute: Route = {
    val route = get {
      pathPrefix(separateOnSlashes(swaggerUiPath)) {
        // when the user hits the doc url, redirect to the index.html with api docs specified on the url
        pathEndOrSingleSlash {
          redirect(
            s"$swaggerUiBaseUrl/$swaggerUiPath/index.html?url=$swaggerUiBaseUrl/$swaggerUiDocsPath",
            StatusCodes.TemporaryRedirect)
        } ~ getFromResourceDirectory(s"META-INF/resources/webjars/swagger-ui/$swaggerUiVersion")
      }
    }
    if (swaggerUiFromRoot) route ~ routeFromRoot else route
  }

}

/**
 * Extends the SwaggerUiHttpService to gets UI configuration values from a provided Typesafe Config.
 */
trait SwaggerUiConfigHttpService extends SwaggerUiHttpService {
  /**
   * @return The swagger UI config.
   */
  def swaggerUiConfig: Config

  override def swaggerUiVersion = swaggerUiConfig.getString("uiVersion")

  abstract override def swaggerUiBaseUrl = swaggerUiConfig.as[Option[String]]("baseUrl").getOrElse(super.swaggerUiBaseUrl)

  abstract override def swaggerUiPath = swaggerUiConfig.as[Option[String]]("uiPath").getOrElse(super.swaggerUiPath)

  abstract override def swaggerUiDocsPath = swaggerUiConfig.as[Option[String]]("docsPath").getOrElse(super.swaggerUiDocsPath)
}

/**
 * An extension of HttpService to serve up a resource containing the swagger api as yaml or json. The resource
 * directory and path on the classpath must match the path for route. The resource can be any file type supported by the
 * swagger UI, but defaults to "yaml". This is an alternative to spray-swagger's SwaggerHttpService.
 */
trait SwaggerResourceHttpService {
  /**
   * @return The directory for the resource under the classpath, and in the url
   */
  def swaggerDirectory: String = "swagger"

  /**
   * @return Name of the service, used to map the documentation resource at "/uiPath/serviceName.resourceType".
   */
  def swaggerServiceName: String

  /**
   * @return The type of the resource, usually "yaml" or "json".
   */
  def swaggerResourceType: String = "yaml"

  /**
   * Swagger UI sends HTTP OPTIONS before ALL requests, and expects a status 200 / OK. When true (the default) the
   * swaggerResourceRoute will return 200 / OK for requests for OPTIONS.
   *
   * See also:
   * - https://github.com/swagger-api/swagger-ui/issues/1209
   * - https://github.com/swagger-api/swagger-ui/issues/161
   * - https://groups.google.com/forum/#!topic/swagger-swaggersocket/S6_I6FBjdZ8
   *
   * @return True if status code 200 should be returned for HTTP OPTIONS requests for the swagger resource.
   */
  def swaggerAllOptionsOk: Boolean = true

  /**
   * @return The path to the swagger docs.
   */
  protected def swaggerDocsPath = s"$swaggerDirectory/$swaggerServiceName.$swaggerResourceType"

  /**
   * @return A route that returns the swagger resource.
   */
  final def swaggerResourceRoute: Route = {
    val swaggerDocsDirective = path(separateOnSlashes(swaggerDocsPath))
    val route = get {
      swaggerDocsDirective {
        // Return /uiPath/serviceName.resourceType from the classpath resources.
        getFromResource(swaggerDocsPath)
      }
    }

    if (swaggerAllOptionsOk) {
      route ~ options {
        // Also return status 200 / OK for all OPTIONS requests.
        complete(StatusCodes.OK)
      }
    } else route
  }
}

/**
 * Extends the SwaggerUiHttpService and SwaggerResourceHttpService to serve up both.
 */
trait SwaggerUiResourceHttpService extends SwaggerUiHttpService with SwaggerResourceHttpService {
  override def swaggerUiDocsPath = swaggerDocsPath

  /**
   * @return A route that redirects to the swagger UI and returns the swagger resource.
   */
  final def swaggerUiResourceRoute: Route = swaggerUiRoute ~ swaggerResourceRoute
}
