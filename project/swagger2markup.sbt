/*
Provides imports used by GenerateRestApiDocs.

Borrowed from https://github.com/21re/sbt-swagger-plugin
Uses latest (at time of coding) s2m 1.3.1
 */
val apacheHttpclientV = "4.3.6"
val asmV = "5.2"
val commonsCodecV = "1.10"
val commonsLang3V = "3.5"
val commonsLoggingV = "1.2"
val guavaV = "18.0"
val jacksonV = "2.8.4"
val javaslangV = "2.0.5"
val plexusUtilsV = "3.0.22"
val slf4jV = "1.7.21"
val swagger2markupV = "1.3.1"

resolvers ++= List(
  "JCenter" at "https://jcenter.bintray.com"
)

libraryDependencies ++= List(
  "io.github.swagger2markup" % "swagger2markup" % swagger2markupV
)

dependencyOverrides ++= List(
  "com.fasterxml.jackson.core" % "jackson-annotations" % jacksonV,
  "com.fasterxml.jackson.core" % "jackson-databind" % jacksonV,
  "com.google.guava" % "guava" % guavaV,
  "commons-codec" % "commons-codec" % commonsCodecV,
  "commons-logging" % "commons-logging" % commonsLoggingV,
  "io.javaslang" % "javaslang" % javaslangV,
  "org.apache.commons" % "commons-lang3" % commonsLang3V,
  "org.apache.httpcomponents" % "httpclient" % apacheHttpclientV,
  "org.codehaus.plexus" % "plexus-utils" % plexusUtilsV,
  "org.ow2.asm" % "asm" % asmV,
  "org.ow2.asm" % "asm-tree" % asmV,
  "org.slf4j" % "slf4j-api" % slf4jV
)
