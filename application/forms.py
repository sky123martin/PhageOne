from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, IntegerField, SubmitField, SelectField
from wtforms.validators import DataRequired, Optional
from flask_wtf.file import FileField, FileRequired, FileAllowed
from flask import request

class BRED_edit_form(FlaskForm):
    phage = StringField ('phage', validators=[DataRequired()])
    # where to insert
    bp_position_lower = IntegerField('start', validators=(Optional(),))
    bp_position_upper = IntegerField('stop', validators=(Optional(),))
    # or
    gp_number = IntegerField('gene product number', validators=(Optional(),))
    # what to swap in
    template_DNA = StringField ('template_DNA', validators=(Optional(),))

    # submit
    submit = SubmitField('Generate BRED Substrates')
