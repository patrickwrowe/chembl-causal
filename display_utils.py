from IPython.display import Image, display


def view_pydot(pydot_graph):
    display(Image(pydot_graph.create_png()))