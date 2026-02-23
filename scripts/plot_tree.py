import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import analysis.event_tree as et
import numpy as np

dict_list = et.read_event_file('nEXO_OD001_fort.72')
df = et.get_dataframe(dict_list)
et.plot_tree_plotly(df, show_geometry= True, show_EM=False)
