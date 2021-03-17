from flask_wtf import FlaskForm
from wtforms import StringField, BooleanField, IntegerField, SubmitField, SelectField
from wtforms.validators import DataRequired, Optional
from flask_wtf.file import FileField, FileRequired, FileAllowed
from flask import request

class BRED_edit_form(FlaskForm):
    phage = StringField ('phage', validators=[DataRequired()])

    # where to insert
    bp_position_start = IntegerField('start', validators=(Optional(),))
    bp_position_stop = IntegerField('stop', validators=(Optional(),))
    # or
    gp_number = IntegerField('gene product number', validators=(Optional(),))

    # what to swap in
    template_DNA = StringField ('template_DNA', validators=(Optional(),))
    # type of edit (replacement, insertion, deletion)
    edit_type = StringField ('edit_type', validators=[DataRequired()])
    location_type = StringField ('location_type', validators=(Optional(),))
    orientation = StringField ('orientation', validators=(Optional(),))
    EGFP = BooleanField('EGFP', validators=(Optional(),))
    # submit
    submit = SubmitField('Generate BRED Substrates')

class editing_guide_form(FlaskForm):
    phage = StringField ('phage', validators=[DataRequired()])

    # TODO allow for editing guides in terms of certain subsets of phages ie clusters or morphotype