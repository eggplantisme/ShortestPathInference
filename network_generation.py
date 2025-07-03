import networkx as nx
from rhg import rhg

def test0():
    G = nx.karate_club_graph()
    path = "./dk-lab-2020_code_hyperlink/netdata/karateclub.net"
    with open(path, 'w') as f:
        for line in nx.generate_edgelist(G, data=False):
            f.write(line+'\n')

def test1():
    n=1000
    R=10
    T=0.5
    alpha=0.75
    hyperbolic_graph = rhg(n, R, T, alpha)
    path = f"./dk-lab-2020_code_hyperlink/netdata/rhg_test_n{n}.net"
    hyperbolic_graph_A = hyperbolic_graph.generate()
    rhg.save(hyperbolic_graph_A, path)


if __name__ == "__main__":
    # test0()
    test1()
