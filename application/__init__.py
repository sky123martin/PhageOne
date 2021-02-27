from flask import Flask
from application.config import Config
import subprocess

app = Flask(__name__)
app.config.from_object(Config)

from flask_bootstrap import Bootstrap

bootstrap = Bootstrap(app)

proc = subprocess.check_call("mkdir -p fasta_files", shell=True)
proc = subprocess.check_call("mkdir -p genes_by_phage", shell=True)


from application import routes