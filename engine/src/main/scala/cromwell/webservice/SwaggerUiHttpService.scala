package cromwell.webservice

import akka.http.scaladsl.model.{HttpResponse, StatusCodes}
import akka.http.scaladsl.server
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server.Route
import akka.stream.scaladsl.Flow
import akka.util.ByteString
import cromwell.webservice.routes.CromwellApiService

/**
 * Serves up the swagger UI from org.webjars/swagger-ui, lightly editing the index.html.
 */
trait SwaggerUiHttpService {

  private lazy val resourceDirectory = s"META-INF/resources/webjars/swagger-ui/${CromwellApiService.swaggerUiVersion}"

  private val serveIndex: server.Route = {
    val swaggerOptions =
      """
        |        validatorUrl: null,
        |        apisSorter: "alpha",
        |        operationsSorter: "alpha"
      """.stripMargin

    mapResponseEntity { entityFromJar =>
      entityFromJar.transformDataBytes(Flow.fromFunction[ByteString, ByteString] { original: ByteString =>
        ByteString(
          original.utf8String
            .replace("""url: "https://petstore.swagger.io/v2/swagger.json"""", "url: 'cromwell.yaml'")
            .replace("""layout: "StandaloneLayout"""", s"""layout: "StandaloneLayout", $swaggerOptions""")
        )
      })
    } {
      getFromResource(s"$resourceDirectory/index.html")
    }
  }

  /**
   * Serves up the swagger UI only. Redirects requests to the root of the UI path to the index.html.
   *
   * @return Route serving the swagger UI.
   */
  final def swaggerUiRoute: Route = {
    pathEndOrSingleSlash {
      get {
        serveIndex
      }
    } ~
      // We have to be explicit about the paths here since we're matching at the root URL and we don't
      // want to catch all paths lest we circumvent Spray's not-found and method-not-allowed error
      // messages.
      (pathPrefixTest("swagger-ui") | pathPrefixTest("oauth2") | pathSuffixTest("js")
        | pathSuffixTest("css") | pathPrefixTest("favicon")) {
        get {
          getFromResourceDirectory(resourceDirectory)
        }
      } ~
      // Redirect legacy `/swagger` or `/swagger/index.html?url=/swagger/cromwell.yaml#fragment` requests to the root
      // URL. The latter form is (somewhat magically) well-behaved in throwing away the `url` query parameter that was
      // the subject of the CVE linked below while preserving any fragment identifiers to scroll to the right spot in
      // the Swagger UI.
      // https://github.com/swagger-api/swagger-ui/security/advisories/GHSA-qrmm-w75w-3wpx
      (path("swagger" / "index.html") | path ("swagger")) {
        get {
          redirect("/", StatusCodes.MovedPermanently)
        }
      }
  }
}
/**
 * An extension of HttpService to serve up a resource containing the swagger api as yaml or json. The resource
 * directory and path on the classpath must match the path for route. The resource can be any file type supported by the
 * swagger UI, but defaults to "yaml". This is an alternative to spray-swagger's SwaggerHttpService.
 */
trait SwaggerResourceHttpService {

  def getBasePathOverride(): Option[String] = {
    Option(System.getenv("SWAGGER_BASE_PATH"))
  }

  /**
   * @return The directory for the resource under the classpath, and in the url
   */
  private val swaggerDirectory: String = "swagger"

  /**
   * @return Name of the service, used to map the documentation resource at "/uiPath/serviceName.resourceType".
   */
  def swaggerServiceName: String

  /**
   * @return The type of the resource, usually "yaml" or "json".
   */
  def swaggerResourceType: String = "yaml"

  /**
   * @return The path to the swagger docs.
   */
  private lazy val swaggerDocsPath = s"$swaggerDirectory/$swaggerServiceName.$swaggerResourceType"

  /**
   * @return A route that returns the swagger resource.
   */
  final def swaggerResourceRoute: Route = {
    // Serve Cromwell API docs from either `/swagger/cromwell.yaml` or just `cromwell.yaml`.
    val swaggerDocsDirective = path(separateOnSlashes(swaggerDocsPath)) | path(s"$swaggerServiceName.$swaggerResourceType")

    def injectBasePath(basePath: Option[String])(response: HttpResponse): HttpResponse = {
      basePath match {
        case _ if response.status != StatusCodes.OK => response
        case None => response
        case Some(base_path) => response.mapEntity { entity =>
          val swapperFlow: Flow[ByteString, ByteString, Any] = Flow[ByteString].map(byteString => ByteString.apply(byteString.utf8String.replace("#basePath: ...", "basePath: " + base_path)))
          entity.transformDataBytes(swapperFlow)
        }
      }
    }

    val route = get {
      swaggerDocsDirective {
        // Return /uiPath/serviceName.resourceType from the classpath resources.
        mapResponse(injectBasePath(getBasePathOverride()))(getFromResource(swaggerDocsPath))
      }
    }

    route ~ options {
      // Also return status 200 / OK for all OPTIONS requests.
      complete(StatusCodes.OK)
    }
  }
}

/**
 * Extends the SwaggerUiHttpService and SwaggerResourceHttpService to serve up both.
 */
trait SwaggerUiResourceHttpService extends SwaggerUiHttpService with SwaggerResourceHttpService {
  /**
   * @return A route that redirects to the swagger UI and returns the swagger resource.
   */
  final def swaggerUiResourceRoute: Route = swaggerUiRoute ~ swaggerResourceRoute
}
