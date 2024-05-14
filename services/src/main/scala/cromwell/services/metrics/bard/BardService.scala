package cromwell.services.metrics.bard

import bio.terra.bard.api.DefaultApi
import bio.terra.bard.client.ApiClient
import bio.terra.bard.model.EventsEventLogRequest
import com.typesafe.scalalogging.LazyLogging
import cromwell.services.metrics.bard.model.BardEvent
import org.apache.http.impl.client.HttpClients
import org.apache.http.impl.conn.PoolingHttpClientConnectionManager
import org.springframework.http.client.HttpComponentsClientHttpRequestFactory
import org.springframework.web.client.RestTemplate

class BardService(bardUrl: String, connectionPoolSize: Int) extends LazyLogging {

  private val restTemplate = makeRestTemplateWithPooling
  private val appId = "cromwell"

  def sendEvent(bearerToken: String, event: BardEvent): Unit = {
    val api = getEventApi(bearerToken)
    val eventLogRequest = new EventsEventLogRequest().properties(event.getProperties)
    try
      api.eventsEventLog(event.eventName, appId, eventLogRequest)
    catch {
      // Sending events to Bard is a best-effort affair. If it fails, log the error and move on.
      case e: Exception =>
        logger.error(s"Failed to send event to Bard: ${e.getMessage}", e)
    }
    ()
  }

  private def getEventApi(bearerToken: String): DefaultApi = {
    val bardClient = new ApiClient(restTemplate)
    bardClient.setBasePath(bardUrl)
    bardClient.setBearerToken(bearerToken)
    new DefaultApi(bardClient)
  }

  /**
    * @return a new RestTemplate backed by a pooling connection manager
    */
  private def makeRestTemplateWithPooling: RestTemplate = {
    val poolingConnManager = new PoolingHttpClientConnectionManager()
    poolingConnManager.setMaxTotal(connectionPoolSize)
    poolingConnManager.setDefaultMaxPerRoute(connectionPoolSize)
    val httpClient = HttpClients.custom.setConnectionManager(poolingConnManager).build
    val factory = new HttpComponentsClientHttpRequestFactory(httpClient)
    new RestTemplate(factory)
  }

}
