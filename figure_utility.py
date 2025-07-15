import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def hyperbolic_graph_visulize(G, r_theta_coords, shortest_path=None, base_node_size=8, shortest_path_node_size=10, shortest_path_edge_size=0.5, base_node_alpha="63", base_edge_alpha="11"):
    if shortest_path is not None:
        src_node = shortest_path[0]
        des_node = shortest_path[-1]
        shortest_path_edges = [(shortest_path[i], shortest_path[i+1]) for i in range(len(shortest_path)-1)]
    pos = dict()
    for i in r_theta_coords.keys():
        r = r_theta_coords[i][0]
        theta = r_theta_coords[i][1]
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        pos[i] = (x, y)
        if shortest_path is not None:
            if i == src_node:
                src_rtheta = (r, theta)
                print(src_rtheta)
            if i == des_node:
                des_rtheta = (r, theta)
    # set node color and size
    node_color = []
    node_size = []
    for node in G.nodes:
        base_node_color = "#91CAE8" + base_node_alpha
        if shortest_path is not None:
            if node == src_node:
                node_color.append("red")
                node_size.append(shortest_path_node_size)
            elif node == des_node:
                node_color.append("green")
                node_size.append(shortest_path_node_size)
            elif node in shortest_path:
                node_color.append("purple")
                node_size.append(shortest_path_node_size)
            else:
                node_color.append(base_node_color)
                node_size.append(base_node_size)
        else:
            node_color.append(base_node_color)
            node_size.append(base_node_size)
    # set edge color and width
    edge_color = []
    width = []
    for edge in G.edges:
        base_edge_color = "#CFDAEC" + base_edge_alpha
        if shortest_path is not None and (edge in shortest_path_edges or (edge[1], edge[0]) in shortest_path_edges):
            edge_color.append("#F3D707")
            width.append(shortest_path_edge_size)
        else:
            edge_color.append(base_edge_color)
            width.append(0.2)
    nx.draw(G, pos=pos, node_size=node_size, node_color=node_color, edge_color=edge_color, width=width)
    plt.axis('equal')
    # visualize geodesic
    if shortest_path is not None:
        KL_src_rtheta = (np.tanh(src_rtheta[0]), src_rtheta[1])
        KL_des_rtheta = (np.tanh(des_rtheta[0]), des_rtheta[1])
        KL_src_xy = (KL_src_rtheta[0]*np.cos(KL_src_rtheta[1]), KL_src_rtheta[0]*np.sin(KL_src_rtheta[1]))
        KL_des_xy = (KL_des_rtheta[0]*np.cos(KL_des_rtheta[1]), KL_des_rtheta[0]*np.sin(KL_des_rtheta[1]))
        k = (KL_src_xy[1]-KL_des_xy[1])/(KL_src_xy[0]-KL_des_xy[0])
        b = KL_des_xy[1] - k * KL_des_xy[0]
        geodesic_x = []
        geodesic_y = []
        for t in np.arange(0, 1.01, 0.01):
            KL_x = (KL_des_xy[0] - KL_src_xy[0]) * t + KL_src_xy[0]
            KL_y = k * KL_x + b
            r = np.arctanh(np.sqrt(KL_x**2+KL_y**2))
            theta = np.arctan2(KL_y, KL_x)  # extend arctan to -pi~pi
            if t == 0:
                print(KL_src_rtheta)
                print(KL_x, KL_y)
                print(r, theta)
            geodesic_x.append(r*np.cos(theta))
            geodesic_y.append(r*np.sin(theta))
        plt.plot(geodesic_x, geodesic_y)
    