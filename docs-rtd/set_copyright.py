from datetime import datetime
def on_config(config, **kwargs):
    config.copyright = f"Copyright © {datetime.now().year} <a href=\"https://broadinstitute.org\">Broad Institute</a>"
