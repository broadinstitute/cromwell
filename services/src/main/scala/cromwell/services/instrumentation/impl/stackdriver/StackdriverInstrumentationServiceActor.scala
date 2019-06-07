package cromwell.services.instrumentation.impl.stackdriver

import java.util

import akka.actor.{Actor, ActorRef, Props}
import com.google.api.gax.core.FixedCredentialsProvider
import com.google.api.{Metric, MetricDescriptor, MonitoredResource}
import com.google.auth.oauth2.GoogleCredentials
import com.google.cloud.monitoring.v3.{MetricServiceClient, MetricServiceSettings}
import com.google.monitoring.v3._
import com.google.protobuf.util.Timestamps
import com.typesafe.config.Config
import cromwell.services.instrumentation.InstrumentationService.InstrumentationServiceMessage
import cromwell.services.instrumentation._
import cromwell.services.instrumentation.impl.stackdriver.StackdriverInstrumentationServiceActor._
import cromwell.services.instrumentation.impl.stackdriver.StackdriverConfig._
import cromwell.util.GracefulShutdownHelper.ShutdownCommand

import scala.collection.JavaConverters._
import scala.concurrent.duration._


class StackdriverInstrumentationServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor {
  implicit lazy val executionContext = context.dispatcher

  final val ActorCreationTime = Timestamps.fromMillis(System.currentTimeMillis())
  final val CustomMetricDomain = "custom.googleapis.com"
  final val InitialDelay = 1.minute

  lazy val stackdriverConfig = StackdriverConfig(serviceConfig)
  lazy val projectName: ProjectName = ProjectName.of(stackdriverConfig.googleProject)
  lazy val metricLabelsMap = generateMetricLabels()

  var metricsMap = Map.empty[StackdriverMetric, List[Double]]

  //TODO: Saloni- add auth in config
  val appDefaultCred = GoogleCredentials.getApplicationDefault

  // Instantiates a client
  val metricServiceSettings = MetricServiceSettings.newBuilder.setCredentialsProvider(FixedCredentialsProvider.create(appDefaultCred)).build
  final val metricServiceClient = MetricServiceClient.create(metricServiceSettings)

  // Prepares the monitored resource descriptor
  lazy val resourceLabels = Map[String, String](("project_id", stackdriverConfig.googleProject))
  lazy val monitoredResource = MonitoredResource.newBuilder.setType("global").putAllLabels(resourceLabels.asJava).build

  // Give the actor time to warm up, then start sending the metrics to Stackdriver at an interval
  context.system.scheduler.schedule(InitialDelay, stackdriverConfig.flushRate, self, SendStackdriverMetricCommand)


  override def receive = {
    case SendStackdriverMetricCommand => sendMetricData()
    case InstrumentationServiceMessage(cromwellMetric) => cromwellMetric match {
      case CromwellTiming(bucket, value, _) => updateMetricMap(bucket, value.toMillis.toDouble, StackdriverGauge)
      case CromwellGauge(bucket, value) => updateMetricMap(bucket, value.toDouble, StackdriverGauge)
      case CromwellCount(bucket, value, _) => updateMetricMap(bucket, value.toDouble, StackdriverCumulative)
      case CromwellIncrement(bucket) => updateMetricMap(bucket, metricValue = 1D, metricKind = StackdriverCumulative)
    }
    case ShutdownCommand =>
      // flush out metrics (if any) before shut down
      sendMetricData()
      metricServiceClient.close()
      context stop self
  }


  private def labelFromConfig(op: StackdriverConfig => Option[String], key: String): List[(String, String)] = {
    op(stackdriverConfig).fold(List.empty[(String, String)])(v => List((key.replace("-", "_"), v)))
  }


  private def generateMetricLabels(): Map[String, String] = {
    val labelsList = labelFromConfig(_.cromwellInstanceIdentifier, CromwellInstanceIdentifier) ++
      labelFromConfig(_.cromwellInstanceRole, CromwellInstanceRole) ++
      labelFromConfig(_.cromwellPerfTestCase, CromwellPerfTest)

    labelsList.toMap
  }


  private def updateMetricMap(bucket: CromwellBucket, metricValue: Double, metricKind: StackdriverMetricKind): Unit = {
    val metricObj = StackdriverMetric(bucket.toStackdriverString, metricKind)

    if (metricsMap.contains(metricObj)) {
      val valueList = metricsMap(metricObj)

      metricsMap += metricObj ->  valueList.::(metricValue)
    }
    else metricsMap += metricObj -> List(metricValue)
  }


  private def sendMetricData(): Unit = {
    metricsMap.foreach { case (key, value: List[Double]) =>
      val dataPointListSum = value.sum

      val dataPoint = key.kind match {
        case StackdriverGauge => dataPointListSum / value.length
        case StackdriverCumulative => dataPointListSum
      }

      writeMetrics(key, dataPoint)
    }

    metricsMap = Map.empty[StackdriverMetric, List[Double]]
  }


  private def timeInterval(metricKind: StackdriverMetricKind): TimeInterval = {
     metricKind match {
       case StackdriverGauge => TimeInterval.newBuilder.setEndTime(Timestamps.fromMillis(System.currentTimeMillis)).build
       case StackdriverCumulative => TimeInterval.newBuilder.setStartTime(ActorCreationTime).setEndTime(Timestamps.fromMillis(System.currentTimeMillis)).build
     }
  }


  private def createTimeSeries(metricKind: StackdriverMetricKind, metric: Metric, resource: MonitoredResource, dataPointList: util.List[Point]): TimeSeries = {
    metricKind match {
      case StackdriverGauge => TimeSeries.newBuilder.setMetric(metric).setResource(resource).addAllPoints(dataPointList).build
      case StackdriverCumulative => TimeSeries.newBuilder.setMetric(metric).setResource(resource).setMetricKind(MetricDescriptor.MetricKind.CUMULATIVE).addAllPoints(dataPointList).build
    }
  }


  private def createMetricBuilder(metricName: String): Metric = {
    metricLabelsMap match {
      case map: Map[String, String] if map.isEmpty => Metric.newBuilder.setType(s"$CustomMetricDomain/$metricName").build
      case map: Map[String, String] if map.nonEmpty => Metric.newBuilder.setType(s"$CustomMetricDomain/$metricName").putAllLabels(map.asJava).build
    }
  }


  private def writeMetrics(metricObj: StackdriverMetric, value: Double): Unit = {
    // Prepares an individual data point
    val interval = timeInterval(metricObj.kind)
    val pointValue = TypedValue.newBuilder().setDoubleValue(value).build()
    val dataPoint: Point = Point.newBuilder.setInterval(interval).setValue(pointValue).build
    val dataPointList: List[Point] = List[Point](dataPoint)

    // Prepares the metric descriptor
    val metric: Metric = createMetricBuilder(metricObj.name)

    // Prepares the time series request
    val timeSeries = createTimeSeries(metricObj.kind, metric, monitoredResource, dataPointList.asJava)
    val timeSeriesList = List[TimeSeries](timeSeries)

    val timeSeriesRequest = CreateTimeSeriesRequest.newBuilder.setName(projectName.toString).addAllTimeSeries(timeSeriesList.asJava).build

    // Writes time series data
    metricServiceClient.createTimeSeries(timeSeriesRequest)

    println(s"Kind: ${metricObj.kind} Value: $value Metric: ${metricObj.name}")
  }
}


object StackdriverInstrumentationServiceActor {
  val CromwellMetricPrefix: String = "cromwell"

  def props(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) = Props(new StackdriverInstrumentationServiceActor(serviceConfig, globalConfig, serviceRegistryActor))

  implicit class CromwellBucketEnhanced(val cromwellBucket: CromwellBucket) extends AnyVal {
    /**
      * Transforms a CromwellBucket to a StatsD path
      */
    def toStackdriverString = (cromwellBucket.prefix ++ cromwellBucket.path.toList).mkString("/").replace(" ", "_")
  }
}


sealed trait StackdriverMetricKind
object StackdriverGauge extends StackdriverMetricKind
object StackdriverCumulative extends StackdriverMetricKind


case class StackdriverMetric(name: String, kind: StackdriverMetricKind)


object SendStackdriverMetricCommand
