Security
========

<!---toc start-->
* [Firecloud](#firecloud)
* [Security by sysadmin](#security)
   * [Multi-tenant](#multi-tenant)
<!---toc end-->

# Firecloud

TODO

# Security by sysadmin
__Warning!__

__This section is community-contributed. It is intended as helpful guidance only, and is not endorsed by the Broad Institute.__

Cromwell running in server mode accepts all connections on the configured webservice port. The simplest way to restrict access is by putting an authenticating proxy server in between users and the cromwell server:
 1. Configure a firewall rule on the cromwell server host to deny access to the webservice port (e.g. 8000) from all addresses except a secure proxy host.
 1. Configure `<YourFavoriteWebProxy>` on the proxy host with `<YourFavoriteAuthMechanism>`, to proxy authenticated traffic from the world to the cromwell server. Using Apache `httpd` web server for example with basic htpassword file-based authentication, the configuration might look something like:

 ```Apache
<Location /cromwell>
    Order deny,allow
    Allow from all
    AuthType Basic
    AuthName "Password Required"
    AuthUserFile /path/to/my/htpasswdfile
    Require user someone someoneelse
    ProxyPass http://101.101.234.567:8000     # address of cromwell server web service
</Location>
```

 1. That's it. Users now hit `http://my.proxy.org/cromwell` with authenticated requests, and they're forwarded to port 8000 on the cromwell server host.

## Multi-tenant
The above scheme extends easily to multiple cromwell instances, for use by different groups within an organization for example. If the instances are running on the same host then each instance should be run as its own dedicated service account user, e.g. `cromwell1`, `cromwell2` etc. so that processes running under one cromwell instance cannot access the files of another; different webservice ports must also be configured. If persistent database storage is being used then each instance should be configured with its own database and database user. The proxy configuration above is extended simply by adding another `Location`:

```Apache
<Location /cromwell1>
    Order deny,allow
    Allow from all
    AuthType Basic
    AuthName "Password Required"
    AuthUserFile /path/to/my/htpasswdfile1
    Require user stillanotherperson andanother
    ProxyPass http://101.101.234.567:8001
</Location>
```

