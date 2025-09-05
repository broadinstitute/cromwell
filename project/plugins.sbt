addSbtPlugin("se.marcuslonnberg" % "sbt-docker" % "1.9.0")
addSbtPlugin("com.eed3si9n" % "sbt-assembly" % "1.1.1")
addSbtPlugin("com.github.sbt" % "sbt-git" % "2.0.0")
addSbtPlugin("org.scoverage" % "sbt-scoverage" % "2.0.4")
addSbtPlugin("com.github.cb372" % "sbt-explicit-dependencies" % "0.2.16")
addSbtPlugin("org.scalameta" % "sbt-scalafmt" % "2.5.2")

// TODO
// This plugin is only used for `sbt publish`, which required GCP auth to push to GAR.
// The plugin itself spews a bunch of errors when added to an environment without GCP auth.
// Only enable the plugin when running `sbt publish`, or only when GCP auth is available.
addSbtPlugin("org.latestbit" % "sbt-gcs-plugin" % "1.14.0")
addDependencyTreePlugin
