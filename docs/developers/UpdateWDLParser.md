This page contains instructions for using [Hermes](https://github.com/scottfrazer/hermes) to regenerate the WDL parser after making changes to the WDL grammar file, which is located at `wdl/transforms/$version/src/main/resources/grammar.hgr`. Note that there is a separate grammar file and parser for each WDL version.

## Install Hermes

Requires Python 3 (tested with Python 3.5). [Pyenv](https://github.com/pyenv/pyenv) is the recommended way to install this older Python version. Newer Python versions may not work correctly.

Hermes can be installed via pip (https://pypi.python.org/pypi/hermes-parser/2.0rc6)

However, It is recommended that Hermes is installed from source into a virtual environment because PyPI isn’t updated very often.

First, clone Hermes and create then activate a new Python virtual environment (in directory ./hermes-venv)
```
$ git clone git@github.com:scottfrazer/hermes.git
$ cd hermes
$ python -m venv hermes-venv
$ source hermes-venv/bin/activate
```

In the Hermes repo, make sure you have the latest ‘develop’ branch checked out. Then, install Hermes into the virtual environment (need to separately run pip installs because the Python 3.5 default index URL doesn't work anymore):
```
$ pip install --index-url https://pypi.org/simple/ moody-templates xtermcolor nose pygments hermes_pygments
$ python setup.py install
$ hermes --help
```

With the virtual environment activated, navigate to `cromwell/wdl/transforms/$version` and run the below command, replacing `$version` with the correct package name (e.g. `biscayne` or `cascades`).:
```
hermes generate src/main/resources/grammar.hgr \
          --language=java \
          --directory=src/main/java \
          --name=wdl \
          --java-package=wdl.$version.parser \
          --java-use-apache-commons --java-imports=org.apache.commons.lang3.StringEscapeUtils \
          --header
```
