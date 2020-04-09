from __future__ import absolute_import
import os
import ee
import math
import networkx as nx
from pprint import pformat
from ee.ee_exception import EEException


class Pipeline:
    def __init__(self,steps,name=None):
        dag = nx.DiGraph(name=name)
        for i,step in enumerate(steps):
            dag.add_node(step.__name__,obj=step,idx=i)
            if callable(step._qa):
                i+1
                dag.add_node('qa',obj,step._qa,idx=i)
            if i > 0:
                dag.add_edge(steps[i-1].__name__,step.__name__)

        self.graph = dag

        return

    def compute(self):
        return

    def show(self):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        pos=nx.spring_layout(self.graph) # positions for all nodes
        nx.draw_networkx_nodes(self.graph,pos,node_size=500)
        nx.draw_networkx_edges(self.graph,pos,width=1,arrows=True,edge_color='k')
        nx.draw_networkx_labels(self.graph,pos,font_size=12,font_family='sans-serif')
        plt.show()
        return
