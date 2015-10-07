addSbtPlugin("com.typesafe.sbt" % "sbt-git" % "0.8.5")

addSbtPlugin("com.github.gseitz" % "sbt-release" % "1.0.0")

addSbtPlugin("org.scoverage" % "sbt-scoverage" % "1.0.4")

// Manually uploaded sbt-coveralls custom build to hybrid ivy/maven, without the Resolver.PluginPattern, due to
// artifactory bug: https://www.jfrog.com/jira/browse/RTFACT-6235
resolvers += Resolver.url("artifactory-plugin-snapshots",
  new URL("https://artifactory.broadinstitute.org/artifactory/plugins-snapshot"))(
    Patterns(true, "[organisation]/[module]/[revision]/[artifact](-[classifier]).[ext]")
  )

// Waiting for https://github.com/scoverage/sbt-coveralls/issues/63 to be released
// Also, when bumping the version, see if this other issue is fixed, and remove the customizations to after_success
// in .travis.yml: https://github.com/scoverage/sbt-coveralls/issues/62
addSbtPlugin("org.scoverage" % "sbt-coveralls" % "1.0.1-c64417d-SNAPSHOT")
