from backpressure_report.lib import log_utils
import pathlib

__log_path = pathlib.Path(__file__).parent.parent / 'resources' / 'alpha_logs.json'
LOG = log_utils.build_log_jsons_from_input_files([__log_path])[0]
