The `server` subcommand on the executable JAR will start an HTTP server which can accept workflow files to run as well as check status and output of existing workflows.

The following sub-sections define which HTTP Requests the web server can accept and what they will return.  Example HTTP requests are given in [HTTPie](https://github.com/jkbrzt/httpie) and [cURL](https://curl.haxx.se/)

### REST API Versions

All web server requests include an API version in the url. The current version is `v1`.