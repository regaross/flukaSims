import event_tree as et
import numpy as np

dict_list = et.read_event_file('../single_runs/nEXO_OD001_fort.72')
df = et.get_dataframe(dict_list)
et.plot_tree_plotly(df)