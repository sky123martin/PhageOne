from flask import Flask
from application.config import Config
import subprocess

app = Flask(__name__)
app.config.from_object(Config)

from flask_bootstrap import Bootstrap

bootstrap = Bootstrap(app)

# print("LETS get")
# if not (os.path.isfile("data/phage_metadata.csv") and os.path.isfile("data/cleaned_gene_list.csv")):
#     print("RIIIIPIN")
#     download_phage_info()

proc = subprocess.check_call("mkdir -p data", shell=True)
proc = subprocess.check_call("mkdir -p data/fasta_files", shell=True)
proc = subprocess.check_call("mkdir -p data/genes_by_phage", shell=True)


from application import routes