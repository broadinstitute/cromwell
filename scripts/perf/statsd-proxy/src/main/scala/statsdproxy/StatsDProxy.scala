package statsdproxy

import java.net.{InetSocketAddress, URI}
import java.nio.file.{Paths, StandardOpenOption}

import cats.effect.IO
import com.typesafe.config.ConfigFactory
import com.typesafe.scalalogging.StrictLogging
import fs2.{Pull, Stream}
import fs2.io.file.pulls.{fromPath, writeAllToFileHandle}
import fs2.io.udp.{AsynchronousSocketGroup, Packet, open}
import net.ceedubs.ficus.Ficus._

/**
  * Proxies incoming udp connection to another udp socket, and write all packets to a file.
  */
object StatsDProxy extends App with StrictLogging {
  implicit val ag = AsynchronousSocketGroup()
  implicit val ec = scala.concurrent.ExecutionContext.global

  private val config = ConfigFactory.load()
  private val proxyHost = config.as[String]("proxy.host")
  private val proxyPort = config.as[Int]("proxy.port")
  private val statsdHost = config.as[String]("statsd.host")
  private val statsdPort = config.as[Int]("statsd.port")
  private val filePath = config.as[String]("output-file")
  
  private val proxyAddress = new InetSocketAddress(proxyHost, proxyPort)
  private val redirectAddress = new InetSocketAddress(statsdHost, statsdPort)

  logger.info(s"Redirecting statsd packets to ${redirectAddress.getHostString}:${redirectAddress.getPort}")
  logger.info(s"Writing statsd packets to $filePath")

  val byteStream: Stream[IO, Byte] = for {
    serverSocket <- open[IO](proxyAddress)
    _ = logger.info(s"Proxy listening at ${proxyAddress.getHostString}:${proxyAddress.getPort}")
    clientSocket <- open[IO]()
    receivedPacket <- serverSocket.reads()
    // Re-wrap the packet in the redirect address, send it, and then unpack it to create a stream of bytes that can be written to a file
    pipedToOutbound <- Stream(Packet(redirectAddress, receivedPacket.bytes))
      .covary[IO]
      .to(clientSocket.writes())
      .map(_ => receivedPacket.bytes.toVector :+ '\n'.toByte)
      .flatMap(Stream.emits(_))
  } yield pipedToOutbound

  val server: Pull[IO, Nothing, Unit] = for {
    filePull <- fromPath[IO](Paths.get(URI.create(filePath)), List(StandardOpenOption.CREATE, StandardOpenOption.WRITE))
    _ <- writeAllToFileHandle(byteStream, filePull.resource)
  } yield ()

  server.stream.compile.toVector.unsafeRunSync()
}
