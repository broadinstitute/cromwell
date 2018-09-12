# File Transfer Protocol (FTP)

Cromwell supports communication with basic FTP servers (not FTPS, and not SFTP either).
Please read the [#Filesystems](Filesystems.md) section before to get an explanation of how filesystems can be configured in general.
For a full example of FTP configuration see [here](#Example)

**Note**: Be aware that FTP is in many cases very inefficient, not secure and generally doesn't scale well. If possible, consider using an alternate file system.

## Overview

Cromwell handles FTP connections as follows:

- Cromwell maintains one FTP FileSystem per FTP server per user
- An FTP FileSystem maintains a pool of connections with a fixed size (see `max-connection-per-server-per-user` configuration below) to that server, for that user
- When an FTP FileSystem hasn't been used in a certain amount of time (see `cache-ttl` configuration below), the associated connections are closed and it is destroyed
- In a given pool, when a connection has been idle for a certain amount of time (see `idle-connection-timeout` configuration below), it is closed

```
 Cromwell
+---------------------------------------------------------------------------------------------------+
|                                                                                                   |
|  FileSystem 1                     FileSystem 2                     FileSystem 3                   |
| +-----------------------------+  +-----------------------------+  +-----------------------------+ |
| |                             |  |                             |  |                             | |
| | ftp://user1@server1.com     |  | ftp://user2@server1.com     |  | ftp://user1@server2.com     | |
| |                             |  |                             |  |                             | |
| |  Connection Pool            |  |  Connection Pool            |  |  Connection Pool            | |
| | +-------------------------+ |  | +-------------------------+ |  | +-------------------------+ | |
| | |                         | |  | |                         | |  | |                         | | |
| | |  - Connection 1         | |  | |  - Connection 1         | |  | |  - Connection 1         | | |
| | |                         | |  | |                         | |  | |                         | | |
| | |  - Connection 2         | |  | |  - Connection 2         | |  | |  - Connection 2         | | |
| | |                         | |  | |                         | |  | |                         | | |
| | |  - Connection 3         | |  | |  - Connection 3         | |  | |  - Connection 3         | | |
| | |                         | |  | |                         | |  | |                         | | |
| | +-------------------------+ |  | +-------------------------+ |  | +-------------------------+ | |
| |                             |  |                             |  |                             | |
| +-----------------------------+  +-----------------------------+  +-----------------------------+ |
+---------------------------------------------------------------------------------------------------+
```


## Global Configuration

The FTP filesystem supports a global configuration. The goal is to configure values that will apply Cromwell-wide, as opposed to other configuration that could be backend (or engine) specific (see below).

Default configuration:

```hocon
filesystems.ftp.global.config = {
    # This value should be sufficiently high to cover for the duration of the longest expected I/O operation.
    # It is a time to live value after which a filesystem (unique per user per server) will be closed and evicted from the cache if unused for the specified duration.
    cache-ttl = 1 day
    
    # How long to wait trying to obtain a connection from the pool before giving up. Don't specify for no timeout
    # obtain-connection-timeout = 1 hour
    
    # Maximum number of connections that will be established per user per ftp server. This is across the entire Cromwell instance.
    # Setting this number allows to workaround FTP server restrictions for the number of connections for a single user per IP address.
    # Has to be >= 2 to allow copying (Copying and FTP file requires downloading and uploading it)
    max-connection-per-server-per-user = 30
    
    # Time after which a connection will be closed if idle. This is to try to free connections from a filesystem when it's not heavily used.
    idle-connection-timeout = 1 hour
    
    # FTP connection port to use
    connection-port: 21
    
    # FTP connection mode
    connection-mode = "passive"
}

```

Note that this configuration being global, it applies the same way regardless of the FTP server. There is currently no good way to make this configuration dependent on the server.

## Instance Configuration

Configuration that can be applied to a backend or engine filesystem stanza:

```hocon
ftp {
  # optional
  auth {
    username = "username"
    password = "password"
    # Optional
    account = "account"
  }
}
```

A default authentication can be provided in the configuration. It can be overridden using the following workflow options: `ftp-username`, `ftp-password`, `ftp-account`
If authentication is omitted the connection is established as an anonymous user.

## Example

Example of configuration mixing the above two sections:

```hocon
filesystems.ftp.global.config {
    # We only override what changes from the default
    max-connection-per-server-per-user = 20
}

# Configure the default auth for the engine
engine.filesystems.ftp {
  auth {
    username = "me"
    password = "my_password"
  }
}

# Configure the default auth for the local backend
backend.providers.Local.config.filesystems.ftp {
  auth {
    username = "me"
    password = "my_password"
  }
}
```
