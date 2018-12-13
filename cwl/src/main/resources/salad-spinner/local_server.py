from flask import Flask
from flask import request
from main import cwltool_salad

app = Flask(__name__)

@app.route('/salad', methods=['GET', 'POST'])
def salad_route():
    return cwltool_salad(request)
