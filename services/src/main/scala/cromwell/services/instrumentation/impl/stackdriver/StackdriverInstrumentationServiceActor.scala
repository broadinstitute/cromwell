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
import cromwell.util.GracefulShutdownHelper.ShutdownCommand
import net.ceedubs.ficus.Ficus._

import scala.collection.JavaConverters._
import scala.concurrent.duration._

sealed trait StackdriverMetricKind

object StackdriverGauge extends StackdriverMetricKind
object StackdriverCumulative extends StackdriverMetricKind


case class StackdriverMetric(name: String, kind: StackdriverMetricKind)

object SendStackdriverMetricCommand


class StackdriverInstrumentationServiceActor(serviceConfig: Config, globalConfig: Config, serviceRegistryActor: ActorRef) extends Actor {
  implicit lazy val executionContext = context.dispatcher

  final val ActorCreationTime = Timestamps.fromMillis(System.currentTimeMillis())
  final val CustomMetricDomain = "custom.googleapis.com"

  var metricsMap = Map.empty[StackdriverMetric, List[Double]]

  //TODO: Saloni- put this in something like StackdriverConfig and please validate them!!
  val projectId = serviceConfig.as[String]("google-project")

  val appDefaultCred = GoogleCredentials.getApplicationDefault

  // Instantiates a client
  val metricServiceSettings = MetricServiceSettings.newBuilder.setCredentialsProvider(FixedCredentialsProvider.create(appDefaultCred)).build
  final val metricServiceClient = MetricServiceClient.create(metricServiceSettings)

  val name = ProjectName.of(projectId)

  context.system.scheduler.schedule(1.minute, 1.minute, self, SendStackdriverMetricCommand)


  override def receive = {
    case SendStackdriverMetricCommand => sendMetricData()
    case InstrumentationServiceMessage(cromwellMetric) => cromwellMetric match {
      case CromwellTiming(bucket, value, _) => updateMetricMap(bucket, value.toMillis.toDouble, StackdriverGauge)
      case CromwellGauge(bucket, value) => updateMetricMap(bucket, value.toDouble, StackdriverGauge)
      case CromwellCount(bucket, value, _) => updateMetricMap(bucket, value.toDouble, StackdriverCumulative)
      case CromwellIncrement(bucket) => updateMetricMap(bucket, metricValue = 1D, metricKind = StackdriverCumulative)
    }
    case ShutdownCommand =>
      metricServiceClient.close()
      context stop self
  }


  private def updateMetricMap(bucket: CromwellBucket, metricValue: Double, metricKind: StackdriverMetricKind) = {
    val metricObj = StackdriverMetric(bucket.toStackdriverString, metricKind)

    if (metricsMap.contains(metricObj)) {
      val valueList = metricsMap(metricObj)

      metricsMap += metricObj ->  valueList.::(metricValue)
    }
    else metricsMap += metricObj -> List(metricValue)

    ()
  }


  private def sendMetricData() = {
    metricsMap.foreach { case (key, value: List[Double]) =>
      val dataPointListSum = value.sum

      val dataPoint = key.kind match {
        case StackdriverGauge => dataPointListSum / value.length
        case StackdriverCumulative => dataPointListSum
      }

      writeMetrics(key, dataPoint)
    }

    println(s"MAP SIZE BEFORE: ${metricsMap.size}")

    metricsMap = Map.empty[StackdriverMetric, List[Double]]

    println(s"MAP SIZE AFTER: ${metricsMap.size}")
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


  private def writeMetrics(metricObj: StackdriverMetric, value: Double) = {
    // Prepares an individual data point
    val interval = timeInterval(metricObj.kind)
    val pointValue = TypedValue.newBuilder().setDoubleValue(value).build()
    val dataPoint: Point = Point.newBuilder.setInterval(interval).setValue(pointValue).build
    val dataPointList: List[Point] = List[Point](dataPoint)

    // Prepares the metric descriptor
    val metric: Metric = Metric.newBuilder.setType(s"$CustomMetricDomain/${metricObj.name}").build

    // Prepares the monitored resource descriptor
    val resourceLabels = Map[String, String](("project_id", projectId))
    val resource: MonitoredResource = MonitoredResource.newBuilder.setType("global").putAllLabels(resourceLabels.asJava).build

    // Prepares the time series request
    val timeSeries = createTimeSeries(metricObj.kind, metric, resource, dataPointList.asJava)
    val timeSeriesList = List[TimeSeries](timeSeries)

    val timeSeriesRequest = CreateTimeSeriesRequest.newBuilder.setName(name.toString).addAllTimeSeries(timeSeriesList.asJava).build

    // Writes time series data
    metricServiceClient.createTimeSeries(timeSeriesRequest)

    println(s"Kind: ${metricObj.kind} Value: $value Metric: ${metricObj.name}")
  }

//  private def writeTimeSeriesMetricsData(bucket: CromwellBucket, metricValue: Double) = {
//    val metricName = (bucket.prefix ++ bucket.path.toList).mkString("/").replace(" ", "_")
//
//    // Prepares an individual data point
//    val interval = TimeInterval.newBuilder.setEndTime(Timestamps.fromMillis(System.currentTimeMillis)).build
//    val pointValue = TypedValue.newBuilder().setDoubleValue(metricValue).build()
//    val dataPoint: Point = Point.newBuilder.setInterval(interval).setValue(pointValue).build
//    val dataPointList = List[Point](dataPoint)
//
//    // Prepares the metric descriptor
//    //    val metricLabels = Map[String, String](("prefix", "Pittsburg"))
//    val metric = Metric.newBuilder.setType(s"$CustomMetricDomain/$metricName").build
//
//    // Prepares the monitored resource descriptor
//    val resourceLabels = Map[String, String](("project_id", projectId))
//    val resource = MonitoredResource.newBuilder.setType("global").putAllLabels(resourceLabels.asJava).build
//
//    // Prepares the time series request
//    val timeSeries = TimeSeries.newBuilder.setMetric(metric).setResource(resource).addAllPoints(dataPointList.asJava).build
//    val timeSeriesList = List[TimeSeries](timeSeries)
//
//    val timeSeriesRequest = CreateTimeSeriesRequest.newBuilder.setName(name.toString).addAllTimeSeries(timeSeriesList.asJava).build
//
//    // Writes time series data
//    metricServiceClient.createTimeSeries(timeSeriesRequest)
//
//    println(s"Time series Metric: $metricName")
//  }
//
//
//  private def writeCumulativeMetricData(bucket: CromwellBucket, count: Long) = {
//    val metricName = (bucket.prefix ++ bucket.path.toList).mkString("/").replace(" ", "_")
//
//    // Prepares an individual data point
//    val interval = TimeInterval.newBuilder.setStartTime(ActorCreationTime).setEndTime(Timestamps.fromMillis(System.currentTimeMillis)).build
//    val pointValue = TypedValue.newBuilder().setDoubleValue(count.toDouble).build()
//    val dataPoint: Point = Point.newBuilder.setInterval(interval).setValue(pointValue).build
//    val dataPointList = List[Point](dataPoint)
//
//    // Prepares the metric descriptor
//    //    val metricLabels = Map[String, String](("prefix", "Pittsburg"))
//    val metric = Metric.newBuilder.setType(s"$CustomMetricDomain/$metricName").build
//
//    // Prepares the monitored resource descriptor
//    val resourceLabels = Map[String, String](("project_id", projectId))
//    val resource = MonitoredResource.newBuilder.setType("global").putAllLabels(resourceLabels.asJava).build
//
//    // Prepares the time series request
//    val timeSeries = TimeSeries.newBuilder.setMetric(metric).setResource(resource).setMetricKind(MetricDescriptor.MetricKind.CUMULATIVE).addAllPoints(dataPointList.asJava).build
//    val timeSeriesList = List[TimeSeries](timeSeries)
//
//    val timeSeriesRequest = CreateTimeSeriesRequest.newBuilder.setName(name.toString).addAllTimeSeries(timeSeriesList.asJava).build
//
//    // Writes time series data
//    metricServiceClient.createTimeSeries(timeSeriesRequest)
//
//    println(s"Cumulative Metric: $metricName")
//  }
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
