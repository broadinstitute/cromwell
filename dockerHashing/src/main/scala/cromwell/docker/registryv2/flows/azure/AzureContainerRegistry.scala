package cromwell.docker.registryv2.flows.azure

import cats.effect.IO
import com.typesafe.scalalogging.LazyLogging
import cromwell.docker.DockerInfoActor.DockerInfoContext
import cromwell.docker.{DockerImageIdentifier, DockerRegistryConfig}
import cromwell.docker.registryv2.DockerRegistryV2Abstract
import org.http4s.{Header, Request}


class AzureContainerRegistry(config: DockerRegistryConfig) extends DockerRegistryV2Abstract(config) with LazyLogging {
  /**
    * (e.g registry-1.docker.io)
    */
  override protected def registryHostName(dockerImageIdentifier: DockerImageIdentifier): String = 
    dockerImageIdentifier.host match {
      case Some(h) => h
      case None => "marketplace.azurecr.io"
    }
  // terradevacrpublic.azurecr.io/terra-azure-relay-listeners:9eb4762
  // mcr.microsoft.com/oss/nginx/nginx:stable

  override def accepts(dockerImageIdentifier: DockerImageIdentifier): Boolean =
    dockerImageIdentifier.hostAsString.contains("azurecr.io") 


  /**
    * (e.g auth.docker.io)
    */
  override protected def authorizationServerHostName(dockerImageIdentifier: DockerImageIdentifier): String =
    dockerImageIdentifier.host match {
      case Some(h) => h
      case None => ???
    }

  override def serviceName: Option[String] =
    throw new Exception("ACR service name is host of user-defined registry, must derive from `DockerImageIdentifier`")
  
  /**
    * Builds the list of headers for the token request
    */
  override protected def buildTokenRequestHeaders(dockerInfoContext: DockerInfoContext): List[Header] = {
    import org.http4s.headers.`Content-Type`
    import org.http4s.MediaType

    List(`Content-Type`(MediaType.application.`x-www-form-urlencoded`))
  }

  /*
  Unlike other repositories, Azure reserves `GET /oauth2/token` for Basic Authentication [0]
  In order to use Oauth we must `POST /oauth2/token` [1]
  
  [0] https://github.com/Azure/acr/blob/main/docs/Token-BasicAuth.md#using-the-token-api
  [1] https://github.com/Azure/acr/blob/main/docs/AAD-OAuth.md#calling-post-oauth2token-to-get-an-acr-access-token 
   */
  override protected def buildTokenRequest(dockerInfoContext: DockerInfoContext): IO[Request[IO]] = {
    import org.http4s.Uri.{Authority, Scheme}
    import org.http4s.client.dsl.io._
    import org.http4s._
    
    val uri = Uri.apply(
      scheme = Option(Scheme.https),
      authority = Option(Authority(host = Uri.RegName(authorizationServerHostName(dockerInfoContext.dockerImageID)))),
      path = "/oauth2/token",
      query = Query.empty
    )

    // val entityBody: EntityBody[F] = EntityEncoder[F, ResourceMetadataRequest].toEntity(body).body
    val request = org.http4s.Method.POST(
      UrlForm(
        "scope" -> "repository:postgres:pull",
        "service" -> dockerInfoContext.dockerImageID.hostAsString,
        "refresh_token" -> "eyJhbGciOiJSUzI1NiIsInR5cCI6IkpXVCIsImtpZCI6IkNPQVU6UERZSDo0SVJYOjM2SEI6TFYzUDpWNFBGOko0NzQ6SzNOSjpPS1JCOlRZQUo6NEc0Szo1Q1NEIn0.eyJqdGkiOiJiNDk4YzJkMS0yMTUwLTQ2NGMtYWQxMy1lN2VkMjNlZjczYzUiLCJzdWIiOiJhbmljaG9sc0BhenVyZS5kZXYuZW52cy10ZXJyYS5iaW8iLCJuYmYiOjE2ODg1OTE3NTQsImV4cCI6MTY4ODYwMzQ1NCwiaWF0IjoxNjg4NTkxNzU0LCJpc3MiOiJBenVyZSBDb250YWluZXIgUmVnaXN0cnkiLCJhdWQiOiJ0ZXJyYWJhdGNoZGV2LmF6dXJlY3IuaW8iLCJ2ZXJzaW9uIjoiMS4wIiwicmlkIjoiYzc0NGM4ZTQwNDNjNDY3N2IxNjQ4NmIxNGU5NTQ4ZDMiLCJncmFudF90eXBlIjoicmVmcmVzaF90b2tlbiIsImFwcGlkIjoiMDRiMDc3OTUtOGRkYi00NjFhLWJiZWUtMDJmOWUxYmY3YjQ2IiwidGVuYW50IjoiZmFkOTA3NTMtMjAyMi00NDU2LTliMGEtYzdlNWI5MzRlNDA4IiwicGVybWlzc2lvbnMiOnsiQWN0aW9ucyI6WyJyZWFkIiwid3JpdGUiLCJkZWxldGUiLCJkZWxldGVkL3JlYWQiLCJkZWxldGVkL3Jlc3RvcmUvYWN0aW9uIl0sIk5vdEFjdGlvbnMiOm51bGx9LCJyb2xlcyI6W119.xRRDFEp3arUC_clpJlgz4D_wqcZK7F9fcDjXMzwbLOFxGkJwW4R4e_lNWL2dDFh-lbBJ95fwywnpbYRyFK3S-csKNCMJe2btACNmX6KVxxM8Ei-bA6AxyBY6OIus94yb3HYJY1a-F5ihA-H4GavyZlkDbMeMc4OIXj90zTrkY-nJHFs0gG8tHWsCfsARRKUsKXpe8Takf8WGyXF9Wy5clhrr3eSlooyl_n4HxJbMa___5tsJG6P_S2A3pMuEYZ4AeB1_4RgPTGxrFc9NFCQ_v2-btdDTAYqhamHKbzlQu7GmvQg9qRUlhYXEcXy8hk9T3pcmzHyluqOcFmyMMMYPgQ",
        "grant_type" -> "refresh_token"
      ),
      uri,
      buildTokenRequestHeaders(dockerInfoContext): _*, // http4s adds `Content-Length` which ACR does not like (400 response) 
    )
    request
  }

  /*
  In Azure, service name does not exist at the registry level, it varies per repo: `terrabatchdev.azurecr.io`
   */
//  override def serviceName(: Option[String] = None
}
