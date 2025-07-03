from shortest_path_finder import ShortestPathFinder
import networkx as nx

p=0.01
t1=0
t2=0
save_path = f"./results/testA1_overlap_perturbate_hgd_results_{1}_{1}_6.26.txt"
ouput_path = f"./dk-lab-2020_code_hyperlink/output/p={p}_{t1}_{t2}_" + save_path.split('/')[-1]
network_path = "./dk-lab-2020_code_hyperlink/netdata/rhg_testA1_n500_alpha0.75_T0.5_d5.net"
perturbate_network_path = f"./dk-lab-2020_code_hyperlink/netdata/perturbate_{p}_" + network_path.split('/')[-1]
edges = []
with open(perturbate_network_path, 'r') as f:
    for line in f.readlines():
        line = line.strip().split(' ')
        edges.append((int(line[0]), int(line[1])))
perturbated_G1 = nx.from_edgelist(edges)
start = 335
end = 273
rtheta_coords = dict()
with open(ouput_path, 'r') as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip().split('\t')
        i = int(line[0])
        r = float(line[2])
        theta = float(line[3])
        rtheta_coords[i] = (r, theta)
length = 4

perturbated_shortest_path = ShortestPathFinder(perturbated_G1).hyperbolic_geodesic_distance(start, end, rtheta_coords=rtheta_coords, length=length)