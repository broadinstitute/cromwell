import json
import logging
import tempfile

from cwltool.loghandler import _logger
from cwltool.load_tool import fetch_document, resolve_tool_uri, validate_document
from cwltool.load_tool import jobloaderctx
from schema_salad.ref_resolver import Loader

def cwltool_salad(request):
    """HTTP Cloud Function.
    Args:
        request (flask.Request): The request object.
        Accepts 2 content types:
        - application/json
            {
              "path": "<http url pointing to a cwl to be saladed>"
            }
        - text/plain or application/octet-stream
            Raw cwl to be saladed
    Returns:
        The saladed version of the cwl passed in the request (or pointed to via url)
    """
    _logger.setLevel(logging.WARN)

    content_type = request.headers['content-type']

    if content_type == 'application/json':
        request_json = request.get_json(silent=True)
        if request_json and 'path' in request_json:
            try:
                saladed = do_salading(request_json['path'])
                return saladed
            except Exception:
                # Fix this to return the error message
                abort(400)
        else:
            raise ValueError("JSON is invalid, or missing a 'path' property")
    elif (content_type == 'application/octet-stream') or (content_type == 'text/plain'):
        temp = tempfile.NamedTemporaryFile()
        try:
            temp.write(request.data)
            temp.seek(0)
            saladed = do_salading(temp.name)
        finally:
            temp.close()
        return saladed
    else:
        raise ValueError("Unknown content type: {}".format(content_type))

def do_salading(path):
    document_loader = Loader(jobloaderctx)
    uri, fileuri = resolve_tool_uri(path, document_loader=document_loader)
    workflowobj = document_loader.fetch(fileuri)
    document_loader, avsc_names, processobj, metadata, uri = validate_document(document_loader, workflowobj, uri, preprocess_only=True)
    return json.dumps(processobj, indent=4)
