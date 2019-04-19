**Warning!**

 - Cromwell is NOT on its own a security appliance!
 - Only YOU are responsible for your own security! 
 - Please be sure to check with your security team before setting up your Cromwell server
 - Some recommendations and suggestions on security can be found below

__This is intended as helpful guidance only, and is not endorsed by the Broad Institute.__

Cromwell running in server mode accepts all connections on the configured webservice port.  Without taking additional measures to protect your Cromwell server, this can leave your Cromwell server and therefore all information stored by and accessible by your Cromwell server vulnerable to anyone who is able to access the server.  For instance, an unprotected Cromwell server running against the [Google Cloud backend](backends/Google) would leave you vulnerable to outside users running workflows that will cost you money, and potentially access the data files you use within your workflows!

The simplest way to restrict access is by putting an authenticating proxy server in between users and the Cromwell server:  

1. Configure a firewall rule on the Cromwell server host to deny access to the webservice port (e.g. 8000) from all addresses except a secure proxy host.  
2. Configure `<YourFavoriteWebProxy>` on the proxy host with `<YourFavoriteAuthMechanism>`, to proxy authenticated traffic from the world to the Cromwell server. Using Apache `httpd` web server for example with basic `htpassword` file-based authentication, the configuration might look something like:

```Apache
<Location /cromwell>
    Order deny,allow
    Allow from all
    AuthType Basic
    AuthName "Password Required"
    AuthUserFile /path/to/my/htpasswdfile
    Require user someone someoneelse
    ProxyPass http://101.101.234.567:8000   # address of cromwell server web service
</Location>
```

Users now hit `http://my.proxy.org/cromwell` with authenticated requests, and they're forwarded to port 8000 on the Cromwell server host. 

**Multiple Servers on one Host**

The above scheme extends easily to multiple Cromwell instances, for use by different groups within an organization, for example. If multiple instances are running on the same host, then each instance should be run as its own dedicated user (such as service account when using a [Google Cloud backend](backends/Google)), e.g. `cromwell1`, `cromwell2` etc so that processes running under one Cromwell instance cannot access the files of another. Different webservice ports must also be configured separately in order to not clash. If persistent database storage is being used, then each instance must be configured with its own database. The proxy configuration above is extended simply by adding another `Location`:

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

**Multiple Tenants in one Cromwell Server**

Even with a proxy in place, a single Cromwell server does not provide _authorization_ for individual users hitting the endpoints.  This allows, for instance, one _authenticated_ user to view the metadata for a workflow run by another _authenticated_ user.  Due to this limitation, it is important that all users of this Cromwell server be trusted.  

With the exception of the use of the `refresh_token` or `google_compute_service_account` scheme in the [Configuration](Configuring#authentication) for the [Google Cloud backend](backends/Google)), Cromwell servers are setup to access files for use in workflows using a single set of credentials.  This means that all users of these Cromwell servers will have access to the same files that any other user has access to.  The `refresh_token` or `google_compute_service_account` scheme is the only way to ensure data is protected among multiple users of a Cromwell server, however the aforementioned caveats of authorization for endpoints still applies.  

**Protecting Secrets**

Various parts of the Cromwell server [Configuration](Configuring) contain sensitive information, e.g. username and password for your persistent database.  It is strongly recommended to protect the configuration files from any untrusted users, for instance by limiting who can access your Cromwell server host or using a technology such as [HashiCorp Vault](https://www.vaultproject.io/).  

Additionally, the contents of the Cromwell database can contain sensitive information, so it is recommended to limit the access of the database only to trusted users.