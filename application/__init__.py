from flask import Flask
from application.config import Config
import subprocess

app = Flask(__name__)
app.config.from_object(Config)

from flask_bootstrap import Bootstrap

bootstrap = Bootstrap(app)

from application import routes