_For the Doc-A-Thon_  
**Questions to answer and things to consider:**

1. Who is visiting the Local page?  
*This is the first in the list of Backends*
2. What do they need to know first?  

3. Is all the important information there? If not, add it!  
*What is an rc file? Write out the full name with the abbreviation, Return Code (rc) file, then abbreviate after.*
4. Are there things that don't need to be there? Remove them.  

5. Are the code and instructions accurate? Try it!

---
 **DELETE ABOVE ONCE COMPLETE**

---


**Local Backend**

The local backend will simply launch a subprocess for each task invocation and wait for it to produce its rc file.

This backend creates three files in the `<call_dir>` (see previous section):

* `script` - A shell script of the job to be run.  This contains the user's command from the `command` section of the WDL code.
* `stdout` - The standard output of the process
* `stderr` - The standard error of the process

The `script` file contains:

```
#!/bin/sh
cd <container_call_root>
<user_command>
echo $? > rc
```

`<container_call_root>` would be equal to `<call_dir>` for non-Docker jobs, or it would be under `/cromwell-executions/<workflow_uuid>/call-<call_name>` if this is running in a Docker container.

When running without docker, the subprocess command that the local backend will launch is:

```
/bin/bash <script>"
```

When running with docker, the subprocess command that the local backend will launch is:

```
docker run --rm -v <cwd>:<docker_cwd> -i <docker_image> /bin/bash < <script>
```

**NOTE**: If you are using the local backend with Docker and Docker Machine on Mac OS X, by default Cromwell can only
run from in any path under your home directory.

The `-v` flag will only work if `<cwd>` is within your home directory because VirtualBox with
Docker Machine only exposes the home directory by default.  Any local path used in `-v` that is not within the user's
home directory will silently be interpreted as references to paths on the VirtualBox VM.  This can manifest in
Cromwell as tasks failing for odd reasons (like missing RC file)

See https://docs.docker.com/engine/userguide/dockervolumes/ for more information on volume mounting in Docker.