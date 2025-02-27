package cromwell.services.metrics.bard

import akka.actor.ActorRef
import bio.terra.bard.api.DefaultApi
import bio.terra.bard.client.ApiClient
import bio.terra.bard.model.EventsEventLogRequest
import cats.data.NonEmptyList
import com.typesafe.scalalogging.LazyLogging
import cromwell.services.instrumentation.CromwellInstrumentation
import cromwell.services.metrics.bard.model.BardEvent
import org.apache.hc.client5.http.impl.classic.HttpClients
import org.apache.hc.client5.http.impl.io.PoolingHttpClientConnectionManager
import org.springframework.http.client.HttpComponentsClientHttpRequestFactory
import org.springframework.web.client.RestTemplate

class BardService(bardUrl: String, connectionPoolSize: Int, serviceRegistry: ActorRef)
    extends LazyLogging
    with CromwellInstrumentation {

  private val restTemplate = makeRestTemplateWithPooling
  private val client = getEventApi(restTemplate)
  private val appId = "cromwell"

  override lazy val serviceRegistryActor: ActorRef = serviceRegistry

  def sendEvent(event: BardEvent): Unit = {
    try {
      val eventLogRequest = new EventsEventLogRequest().properties(event.getProperties)
      client.eventsEventLog(event.eventName, appId, eventLogRequest)
      increment(NonEmptyList.of("send_event", "success"), Some("bard"))
    } catch {
      case e: Exception =>
        logger.error(s"Failed to send event to Bard: ${e.getMessage}", e)
        increment(NonEmptyList.of("send_event", "failure"), Some("bard"))
    }
    ()
  }

  private def getEventApi(restTemplate: RestTemplate): DefaultApi = {
    val bardClient = new ApiClient(restTemplate)
    bardClient.setBasePath(bardUrl)
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
