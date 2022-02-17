package cromiam.webservice

import akka.http.scaladsl.model._
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server._
import akka.stream.Materializer
import akka.stream.scaladsl.Sink
import akka.util.ByteString
import cromiam.server.config.SwaggerOauthConfig

import scala.concurrent.{ExecutionContext, Future}

/**
 * Serves up the swagger UI from org.webjars/swagger-ui.
  *
  * This code is NOT original, and is a copy-paste frankenstein of cromwell, rawls, etc until we create a library.
  *
  * https://github.com/broadinstitute/cromwell/blob/644967de35ad982894b003dcd9ccd6441b76d258/engine/src/main/scala/cromwell/webservice/SwaggerUiHttpService.scala
  * https://github.com/broadinstitute/rawls/blob/4a47a017e0e2a0489a75ef80e02f8f4ff003734a/core/src/main/scala/org/broadinstitute/dsde/rawls/webservice/RawlsApiService.scala#L85-L96
  * https://vimeo.com/215325495
 */
trait SwaggerUiHttpService extends Directives {
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
  def swaggerUiBaseUrl = ""

  /**
   * @return The path to the swagger UI html documents. For example "swagger"
   */
  def swaggerUiPath = "swagger"

  /**
   * The path to the actual swagger documentation in either yaml or json, to be rendered by the swagger UI html.
   *
   * @return The path to the api documentation to render in the swagger UI.
   *         For example "api-docs" or "swagger/lenthall.yaml".
   */
  def swaggerUiDocsPath = "api-docs"

  /**
   * @return When true, if someone requests / (or /baseUrl if setup), redirect to the swagger UI.
   */
  def swaggerUiFromRoot = true

  private def routeFromRoot = get {
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
  final def swaggerUiRoute = {
    // when the user hits the doc url, redirect to the index.html with api docs specified on the url
    val indexRedirect: Route = pathEndOrSingleSlash {
      redirect(
        s"$swaggerUiBaseUrl/$swaggerUiPath/index.html?url=$swaggerUiBaseUrl/$swaggerUiDocsPath",
        StatusCodes.TemporaryRedirect)
    }

    /* Serve a resource from the swagger-ui webjar/bundle */
    val resourceServe: Route = getFromResourceDirectory(s"META-INF/resources/webjars/swagger-ui/$swaggerUiVersion")

    /* Mashup of mapResponseWith and mapResponseEntity */
    def mapResponseEntityWith(f: ResponseEntity => Future[ResponseEntity]): Directive0 = {
      extractExecutionContext flatMap { implicit executionContext =>
        mapRouteResultWithPF {
          case RouteResult.Complete(response) =>
            f(response.entity) map { entity =>
              RouteResult.Complete(response.withEntity(entity))
            }
        }
      }
    }

    /* Serve up the index.html, after passing it through a function that rewrites the response. */
    val indexServe: Route = {
      pathPrefixTest("index.html") {
        extractExecutionContext { implicit executionContext =>
          extractMaterializer { implicit materializer =>
            mapResponseEntityWith(rewriteSwaggerIndex) {
              resourceServe
            }
          }
        }
      }
    }

    /* Redirect to the index, serve the rewritten index, or a resource from the swaggerUI webjar. */
    val route = get {
      pathPrefix(separateOnSlashes(swaggerUiPath)) {
        concat(indexRedirect, indexServe, resourceServe)
      }
    }
    if (swaggerUiFromRoot) route ~ routeFromRoot else route
  }

  private[this] final def rewriteSwaggerIndex(responseEntity: ResponseEntity)
                                             (implicit
                                              executionContext: ExecutionContext,
                                              materializer: Materializer
                                             ): Future[ResponseEntity] = {
    val contentType = responseEntity.contentType
    for {
      data <- responseEntity.dataBytes.runWith(Sink.head) // Similar to responseEntity.toStrict, but without a timeout
      replaced = rewriteSwaggerIndex(data.utf8String)
    } yield HttpEntity.Strict(contentType, ByteString(replaced))
  }

  /** Rewrite the swagger index.html. Default passes through the origin data. */
  protected def rewriteSwaggerIndex(data: String): String = data
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
  def swaggerDirectory = "swagger"

  /**
   * @return Name of the service, used to map the documentation resource at "/uiPath/serviceName.resourceType".
   */
  def swaggerServiceName: String

  /**
   * @return The type of the resource, usually "yaml" or "json".
   */
  def swaggerResourceType = "yaml"

  /**
   * @return The path to the swagger docs.
   */
  protected def swaggerDocsPath = s"$swaggerDirectory/$swaggerServiceName.$swaggerResourceType"

  /**
   * @return A route that returns the swagger resource.
   */
  final def swaggerResourceRoute = {
    // Serve CromIAM API docs from either `/swagger/cromiam.yaml` or just `cromiam.yaml`.
    val swaggerDocsDirective = path(separateOnSlashes(swaggerDocsPath)) | path(s"$swaggerServiceName.$swaggerResourceType")
    val route = get {
      swaggerDocsDirective {
        // Return /uiPath/serviceName.resourceType from the classpath resources.
        getFromResource(swaggerDocsPath)
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
  override def swaggerUiDocsPath = swaggerDocsPath

  def oauthConfig: SwaggerOauthConfig

  /**
   * @return A route that redirects to the swagger UI and returns the swagger resource.
   */
  final def swaggerUiResourceRoute = swaggerUiRoute ~ swaggerResourceRoute

  override protected def rewriteSwaggerIndex(data: String): String = {
    // via https://github.com/swagger-api/swagger-ui/tree/v3.2.2#swaggeruibundle

    val bundleOriginal       = """url: "http://petstore.swagger.io/v2/swagger.json","""
    val bundleOriginalSecure = """url: "https://petstore.swagger.io/v2/swagger.json","""

    val bundleReplacement =
      s"""|url: "/$swaggerDocsPath",
          |validatorUrl: null,
          |oauth2RedirectUrl: window.location.origin + "/$swaggerDirectory/oauth2-redirect.html",
          |tagsSorter: "alpha",
          |operationsSorter: "alpha",
          |""".stripMargin

    // NOTE: Until engine and cromiam can re-common-ize the swagger code, this needs to be sync'ed with
    //   cromiam.webservice.BasicSwaggerUiHttpServiceSpec.rewriteSwaggerIndex
    val initOAuthOriginal = "window.ui = ui"

    val initOAuthReplacement =
      s"""|
          |ui.initOAuth({
          |    clientId: "${oauthConfig.clientId}",
          |    realm: "${oauthConfig.realm}",
          |    appName: "${oauthConfig.appName}",
          |    scopeSeparator: " "
          |  })
          |
          |$initOAuthOriginal
          |""".stripMargin

    data
      .replace(initOAuthOriginal, initOAuthReplacement)
      .replace(bundleOriginal, bundleReplacement)
      .replace(bundleOriginalSecure, bundleReplacement)
  }
}
