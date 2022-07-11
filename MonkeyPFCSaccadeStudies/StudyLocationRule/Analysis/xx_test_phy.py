import os
from phy.apps.template import template_gui

preproc_dir = os.path.join(os.path.dirname(__file__), '..', 'Data', 'Preprocessed')
sess_name = 'sra3_1_j_051_00+'
phy_dir = os.path.join(preproc_dir, sess_name, 'phy')
os.chdir(phy_dir)

template_gui("params.py")
