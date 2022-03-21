package cromiam.webservice

import akka.http.scaladsl.model._
import akka.http.scaladsl.server
import akka.http.scaladsl.server.Directives._
import akka.http.scaladsl.server._
import akka.stream.scaladsl.Flow
import akka.util.ByteString
import common.util.VersionUtil
import cromiam.server.config.SwaggerOauthConfig

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

  def oauthConfig: SwaggerOauthConfig

  private lazy val resourceDirectory = {
    val swaggerUiVersion = VersionUtil.getVersion("swagger-ui", VersionUtil.sbtDependencyVersion("swaggerUi"))
    s"META-INF/resources/webjars/swagger-ui/$swaggerUiVersion"
  }

  private val serveIndex: server.Route = {
    mapResponseEntity { entityFromJar =>
      entityFromJar.transformDataBytes(Flow.fromFunction[ByteString, ByteString] { original: ByteString =>
        ByteString(rewriteSwaggerIndex(original.utf8String))
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
      (path("swagger" / "index.html") | path("swagger")) {
        get {
          redirect("/", StatusCodes.MovedPermanently)

        }
      }
  }

  /** Rewrite the swagger index.html. Default passes through the origin data. */
  protected def rewriteSwaggerIndex(original: String): String = {
    val swaggerOptions =
      s"""
        |        validatorUrl: null,
        |        apisSorter: "alpha",
        |        oauth2RedirectUrl: window.location.origin + "/swagger/oauth2-redirect.html",
        |        operationsSorter: "alpha"
      """.stripMargin

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


    original
      .replace(initOAuthOriginal, initOAuthReplacement)
      .replace("""url: "https://petstore.swagger.io/v2/swagger.json"""", "url: 'cromiam.yaml'")
      .replace("""layout: "StandaloneLayout"""", s"""layout: "StandaloneLayout", $swaggerOptions""")

  }
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
  private lazy val swaggerDocsPath = s"$swaggerDirectory/$swaggerServiceName.$swaggerResourceType"

  /**
   * @return A route that returns the swagger resource.
   */
  final def swaggerResourceRoute: Route = {
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

  /**
   * @return A route that redirects to the swagger UI and returns the swagger resource.
   */
  final def swaggerUiResourceRoute: Route = swaggerUiRoute ~ swaggerResourceRoute
}
