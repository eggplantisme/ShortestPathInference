import networkx as nx
from shortest_path_finder import *
from rhg import *
import random
from multiprocessing import Pool
import os

def read(path, add_paths=None, method="dijkstra", metric="overlap"):
    with open(path, 'r') as f:
        results_list = [row.strip().split() for row in f.readlines()]
        # print(results_list)
        if add_paths is not None:
            print(f"Additional result in {add_paths}")
            for add_path in add_paths:
                with open(add_path, 'r') as add_f:
                    results_list = results_list + [row.strip().split() for row in add_f.readlines()]
        results = np.round(np.float64(results_list), decimals=5)
        # print(results)
        ps = np.unique(results[:, 0])
        overlaps = None
        edit_distance = None
        if metric == "all" or metric == "overlap":
            if method == "all":
                overlaps = np.zeros((np.size(ps), 3))
            else:
                overlaps = np.zeros((np.size(ps), 1))
        if metric == "all" or metric =="edit_distance":
            if method == "all":
                edit_distance = np.zeros((np.size(ps), 3))
            else:
                edit_distance = np.zeros((np.size(ps), 1))
        i = 0
        for p in ps:
            filter_results = results[np.squeeze(np.argwhere(results[:, 0] == p))]
            if np.size(filter_results) == 0:
                print(f"Some parameter p={p} didn't run!")
                i += 1
                continue
            if np.all(filter_results[:, 3:]==-1):
                print(f"p={p} all perturbating disconnect st!")
                if overlaps is not None:
                    overlaps[i, :] = -1
                if edit_distance is not None:
                    edit_distance[i, :] = -1
                i += 1
                continue
            # print(filter_results)
            filter_results = filter_results[np.squeeze(np.argwhere(filter_results[:, 3] != -1))]
            # print(filter_results)
            if len(np.shape(filter_results)) == 1:
                mean_filter_results = filter_results[3:]
            else:
                mean_filter_results = np.mean(filter_results, 0)[3:]
            # print(mean_filter_results)
            if overlaps is not None:
                if method == "all" and metric == "all":
                    overlaps[i, :] = [mean_filter_results[0], mean_filter_results[2], mean_filter_results[4]]
                elif metric == 'all':
                    overlaps[i, :] = mean_filter_results[0]
                    pass
            if edit_distance is not None:
                if method == "all" and metric == "all":
                    edit_distance[i, :] = [mean_filter_results[1], mean_filter_results[3], mean_filter_results[5]]
                elif metric == 'all':
                    edit_distance[i, :] = mean_filter_results[1]
                    pass
            # overlaps[i] = mean_filter_results[3]
            i += 1
    return ps, overlaps, edit_distance

def print_error(value):
    print(value)

def write_results(arg):
    """
    :param arg: savepath, results
    :return:
    """
    if arg[0] is not None:
        with open(arg[0], 'a') as fw:
            fw.write(arg[1])

def calc_metric(metric, perturbated_shortest_path, actual_shortest_paths):
    result = ""
    if metric == "all" or metric == "overlap":
        overlap = ShortestPathFinder.overlap(perturbated_shortest_path, actual_shortest_paths)
        result += f" {overlap}"
    if metric == "all" or metric == "edit_distance":
        edit_distance = ShortestPathFinder.path_edit_similarity(perturbated_shortest_path, actual_shortest_paths)
        result += f" {edit_distance}"
    return result

def run(G, p, shortest_path_times, perturbation_times, save_path, method=None, rtheta_coords=None, network_path=None, gamma=None, T=None, R=None, metric="overlap"):
    # G is fully connected
    result = ""
    edges = G.edges
    nodes = G.nodes
    for t1 in range(shortest_path_times):
        # sample shortest path
        temp = random.sample(list(G.nodes), 2)
        start, end = temp[0], temp[1]
        actual_shortest_paths = [path for path in nx.all_shortest_paths(G, source=start, target=end)]
        print(f"Perturbation p={p}, {t1}th random select {len(actual_shortest_paths)} shortest path, between {start} and {end}, with length {len(actual_shortest_paths[0])}")
        for t2 in tqdm(range(perturbation_times), desc="Perturbating..."):
            """ make sure the start and end in thie GCC otherwise re perturbate (to a maximum times of re perturbation) """
            result += f"{p} {t1} {t2}"
            gcc = set()
            max_try = 10
            i = 0
            while i <= max_try and (start not in gcc or end not in gcc):
                perturbated_edges = rhg.random_perturbation_byedges(edges, p)
                perturbated_G1 = nx.from_edgelist(perturbated_edges)
                # perturbated_G1.add_nodes_from(nodes)  # ensure all nodes are present
                gcc = max(nx.connected_components(perturbated_G1), key=len)
                i += 1
            if i > max_try and (start not in gcc or end not in gcc):
                print(f"p={p} can not keep both start and end in giant connect component!")
                perturbated_shortest_path = None
                if method == "all":
                    # add three method null result
                    for j in range(3):
                        result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
                else:
                    result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
            else:
                try:
                    rtheta_coords = None
                    gcc_perturbated_G1 = perturbated_G1.subgraph(gcc).copy()
                    if method == "all" or method is None or method == "dijkstra":
                        perturbated_shortest_path = ShortestPathFinder(perturbated_G1).dijkstra(start, end)
                        result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
                    if method == "hyperbolic_geodesic_distance":
                        length = len(actual_shortest_paths[0])
                        perturbated_shortest_path = ShortestPathFinder(perturbated_G1).hyperbolic_geodesic_distance(start, end, rtheta_coords=rtheta_coords, length=length)
                        result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
                    if method == "all" or method == "perturbated_hyperbolic_geodesic_distance":
                        length = len(actual_shortest_paths[0])
                        """ 
                            note the node mapping!!!
                        """
                        nodes = perturbated_G1.nodes
                        ouput_path = f"./dk-lab-2020_code_hyperlink/output/p={p}_{t1}_{t2}_" + save_path.split('/')[-1]
                        # if os.path.exists(ouput_path):
                        #     print(f"load exist embedding coords.")
                        #     # print(G.nodes)
                        #     rtheta_coords = rhg.load_hyperlink_coords(ouput_path)
                        # else:
                        perturbate_network_path = f"./dk-lab-2020_code_hyperlink/netdata/perturbate_{p}_" + network_path.split('/')[-1]
                        rhg.save_byedges(perturbated_edges, perturbate_network_path)
                        if rtheta_coords is None:
                            rtheta_coords = rhg.hyperlink_embedding(perturbate_network_path, gamma, T, ouput_path, p=p)
                        perturbated_shortest_path = ShortestPathFinder(gcc_perturbated_G1).hyperbolic_geodesic_distance(start, end, rtheta_coords=rtheta_coords, length=length)
                        result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
                    if method == "all" or method == "perturbated_hyperbolic_p_dijkstra":
                        ouput_path = f"./dk-lab-2020_code_hyperlink/output/p={p}_{t1}_{t2}_" + save_path.split('/')[-1]
                        perturbate_network_path = f"./dk-lab-2020_code_hyperlink/netdata/perturbate_{p}_" + network_path.split('/')[-1]
                        rhg.save_byedges(perturbated_edges, perturbate_network_path)
                        if rtheta_coords is None:
                            rtheta_coords = rhg.hyperlink_embedding(perturbate_network_path, gamma, T, ouput_path, p=p)
                        perturbated_shortest_path = ShortestPathFinder(gcc_perturbated_G1).dijkstra_weighted_rhg_p(start, end, rtheta_coords, R, T)
                        result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
                except nx.exception.NetworkXNoPath:
                    print(f"\nno path between {start}, {end} after remove {p} propotion edges!")
            result += "\n"
            print(result)
            # if perturbated_shortest_path is None:
            #     overlap = 0
            #     edit_path = 0
            #     result += f"{p} {t1} {t2} {overlap}\n"
            # else:
            #     if metric == "overlap":
            #         overlap = ShortestPathFinder.overlap(perturbated_shortest_path, actual_shortest_paths)
            #         result += f"{p} {t1} {t2} {overlap}\n"
            #     elif metric == "edit_distance":
            #         edit_distance = ShortestPathFinder.path_edit_similarity(perturbated_shortest_path, actual_shortest_paths)
            #         result += f"{p} {t1} {t2} {edit_distance}\n"
            #         print("\nedit_distance:", perturbated_shortest_path, actual_shortest_paths, edit_distance)
    return save_path, result

def run_1(G, p, shortest_path_times, perturbation_times, save_path, method=None, given_rtheta_coords=None, network_path=None, gamma=None, T=None, R=None, metric="overlap"):
    # first perturbate then random select shortest path
    # G is fully connected
    result = ""
    edges = G.edges
    nodes = G.nodes
    for t2 in tqdm(range(perturbation_times), desc="Perturbating..."):
        perturbated_edges = rhg.random_perturbation_byedges(edges, p)
        perturbated_G1 = nx.from_edgelist(perturbated_edges)
        # perturbated_G1.add_nodes_from(nodes)  # ensure all nodes are present
        gcc = max(nx.connected_components(perturbated_G1), key=len)
        gcc_perturbated_G1 = perturbated_G1.subgraph(gcc).copy()
        rtheta_coords = None
        mercator_rtheta_coords = None
        for t1 in range(shortest_path_times):
            result += f"{p} {t1} {t2}"
            temp = random.sample(list(gcc_perturbated_G1.nodes), 2)
            start, end = temp[0], temp[1]
            actual_shortest_paths = [path for path in nx.all_shortest_paths(G, source=start, target=end)]
            print(f"p={p}, {t2}th perturbation, {t1}th random select {len(actual_shortest_paths)} shortest path, between {start} and {end}, with length {len(actual_shortest_paths[0])}")
            try:
                if "all" in method or method is None or "dijkstra" in method:
                    perturbated_shortest_path = ShortestPathFinder(perturbated_G1).dijkstra(start, end)
                    result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
                if "hyperbolic_geodesic_distance" in method:
                    length = len(actual_shortest_paths[0])
                    perturbated_shortest_path = ShortestPathFinder(perturbated_G1).hyperbolic_geodesic_distance(start, end, rtheta_coords=given_rtheta_coords, length=length)
                    result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
                if "all" in method or "perturbated_mercator_geodesic_distance" in method:
                    length = len(actual_shortest_paths[0])
                    ouput_path = f"./mercator-master/output/p={p}_{t1}_{t2}_" + save_path.split('/')[-1]
                    perturbate_network_path = f"./mercator-master/netdata/perturbate_{p}_" + network_path.split('/')[-1]
                    gcc_perturbated_edges = nx.to_edgelist(gcc_perturbated_G1)
                    rhg.save_byedges(gcc_perturbated_edges, perturbate_network_path)
                    if mercator_rtheta_coords is None:
                        mercator_rtheta_coords = rhg.mercator_embedding(perturbate_network_path, ouput_path, remove_output_file=True)
                    perturbated_shortest_path = ShortestPathFinder(gcc_perturbated_G1).hyperbolic_geodesic_distance(start, end, rtheta_coords=mercator_rtheta_coords, length=length)
                    result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
                if "all" in method or "perturbated_hyperbolic_geodesic_distance" in method:
                    length = len(actual_shortest_paths[0])
                    """ 
                        note the node mapping!!!
                    """
                    nodes = perturbated_G1.nodes
                    ouput_path = f"./dk-lab-2020_code_hyperlink/output/p={p}_{t1}_{t2}_" + save_path.split('/')[-1]
                    # if os.path.exists(ouput_path):
                    #     print(f"load exist embedding coords.")
                    #     # print(G.nodes)
                    #     rtheta_coords = rhg.load_hyperlink_coords(ouput_path)
                    # else:
                    perturbate_network_path = f"./dk-lab-2020_code_hyperlink/netdata/perturbate_{p}_" + network_path.split('/')[-1]
                    rhg.save_byedges(perturbated_edges, perturbate_network_path)
                    if rtheta_coords is None:
                        rtheta_coords = rhg.hyperlink_embedding(perturbate_network_path, gamma, T, ouput_path, p=p)
                    perturbated_shortest_path = ShortestPathFinder(gcc_perturbated_G1).hyperbolic_geodesic_distance(start, end, rtheta_coords=rtheta_coords, length=length)
                    result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
                if "all" in method or "perturbated_hyperbolic_p_dijkstra" in method:
                    ouput_path = f"./dk-lab-2020_code_hyperlink/output/p={p}_{t1}_{t2}_" + save_path.split('/')[-1]
                    perturbate_network_path = f"./dk-lab-2020_code_hyperlink/netdata/perturbate_{p}_" + network_path.split('/')[-1]
                    rhg.save_byedges(perturbated_edges, perturbate_network_path)
                    if rtheta_coords is None:
                        rtheta_coords = rhg.hyperlink_embedding(perturbate_network_path, gamma, T, ouput_path, p=p)
                    perturbated_shortest_path = ShortestPathFinder(gcc_perturbated_G1).dijkstra_weighted_rhg_p(start, end, rtheta_coords, R, T)
                    result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
                if "all" in method or "perturbated_mercator_p_dijkstra" in method:
                    ouput_path = f"./mercator-master/output/p={p}_{t1}_{t2}_" + save_path.split('/')[-1]
                    perturbate_network_path = f"./mercator-master/netdata/perturbate_{p}_" + network_path.split('/')[-1]
                    gcc_perturbated_edges = nx.to_edgelist(gcc_perturbated_G1)
                    rhg.save_byedges(gcc_perturbated_edges, perturbate_network_path)
                    if mercator_rtheta_coords is None:
                        mercator_rtheta_coords = rhg.mercator_embedding(perturbate_network_path, ouput_path, remove_output_file=True)
                    perturbated_shortest_path = ShortestPathFinder(gcc_perturbated_G1).dijkstra_weighted_rhg_p(start, end, rtheta_coords, R, T)
                    result += calc_metric(metric, perturbated_shortest_path, actual_shortest_paths)
            except nx.exception.NetworkXNoPath:
                print(f"\nno path between {start}, {end} after remove {p} propotion edges!")
            result += "\n"
        # print(result)
    return save_path, result


class exp:
    def __init__(self):
        self.net_path = None
        self.G = None
        self.max_component_subgraph = None
        self.precision_p = 5
        self.ps = None
        self.shortest_path_times = 2
        self.perturbation_times = 2
        self.metric = "overlap"
        self.save_path = None
        self.multiprocessing = False
        self.method = None
        self.rtheta_coords = None
        self.gamma = None
        self.T = None
        self.R = None
    
    def load_network(self, net_path):
        # net1_path = "./dk-lab-2020_code_hyperlink/netdata/h2v1.dat"
        self.net_path = net_path
        self.G = rhg.load(net_path)
        print(f"Network has nodes {len(self.G.nodes)}, edges {len(self.G.edges)}")
        largest_cc = max(nx.connected_components(self.G), key=len)
        print(f"Largest component size {len(largest_cc)}")
        self.max_component_subgraph = self.G.subgraph(largest_cc).copy()
    
    def save_run(self):
        """ find done p """
        p_done = set()
        if os.path.exists(self.save_path):
            with open(self.save_path, 'r') as f:
                for row in f.readlines():
                    row = row.strip().split()
                    p_done.add(np.around(float(row[0]), self.precision_p))
        if self.multiprocessing:
            pool = Pool(2)
            for i, p in enumerate(self.ps):
                if round(p, self.precision_p) in p_done:
                    print(f'Perturbation {p} has been run!')
                    continue
                pool.apply_async(run_1, args=(self.max_component_subgraph, p, self.shortest_path_times, self.perturbation_times, self.save_path, self.method,
                self.rtheta_coords, self.net_path, self.gamma, self.T, self.R, self.metric, ), 
                callback=write_results, error_callback=print_error)
            pool.close()
            pool.join()
        else:
            for i, p in enumerate(self.ps):
                result = ""
                if round(p, self.precision_p) in p_done:
                    print(f'Perturbation {p} has been run!')
                    continue
                _, text_result = run_1(self.max_component_subgraph, p, self.shortest_path_times, self.perturbation_times, self.save_path, self.method, 
                self.rtheta_coords, self.net_path, self.gamma, self.T, self.R, self.metric)
                result += text_result
                write_results((self.save_path, result))


def main0():
    """ prepare network """
    net1_path = "./dk-lab-2020_code_hyperlink/netdata/h2v1.dat"
    G1 = rhg.load(net1_path)
    print(f"Network 1 has nodes {len(G1.nodes)}, edges {len(G1.edges)}")
    largest_cc = max(nx.connected_components(G1), key=len)
    print(f"Largest component size {len(largest_cc)}")
    max_component_subgraph = G1.subgraph(largest_cc).copy()
    """ prepare exp parameter """
    precision_p = 5
    ps = np.around(np.arange(0.1, 1, 0.001), precision_p)
    # ps = np.concatenate((np.linspace(0.0001, 0.1, 10), np.linspace(0.5, 0.9, 2)), axis=None)
    print(f"Number of p {np.size(ps)}")
    shortest_path_times = 2  # 5
    perturbation_times = 2  # 5
    # save_path = "./results/network1_overlap_results_6.24.txt"
    # save_path = "./results/network1_overlap_hgd_results_6.25.txt"
    # save_path = "./results/network1_overlap_dijkstra_results_6.25.txt"
    # save_path = f"./results/network1_overlap_dijkstra_results_{shortest_path_times}_{perturbation_times}_6.26.txt"
    save_path = f"./results/network1_overlap_perturbate_hgd_results_{shortest_path_times}_{perturbation_times}_6.26.txt"

    multiprocessing = False
    # method = "hyperbolic_geodesic_distance"
    # method = "dijkstra"
    method = "perturbated_hyperbolic_geodesic_distance"
    rtheta_coords = None
    gamma = None
    T = None
    if method == "hyperbolic_geodesic_distance":
        # use given actual hyperbolic coordinates
        print(f"For method {method}, get r_theta coordinates.")
        rtheta_coords = dict()
        with open(f"./dk-lab-2020_code_hyperlink/coords_real/h2v1.coords.dat", 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                i = int(line[0])
                r = float(line[1])
                theta = float(line[2])
                rtheta_coords[i] = (r, theta)
    elif method == "perturbated_hyperbolic_geodesic_distance":
        # use hyperlink to embed hyperbolic coordinates in perturbated network, given gamma, T and p
        gamma = 2.5
        T = 0.1
    """ find done p """
    p_done = set()
    if os.path.exists(save_path):
        with open(save_path, 'r') as f:
            for row in f.readlines():
                row = row.strip().split()
                p_done.add(np.around(float(row[0]), precision_p))
    if multiprocessing:
        pool = Pool(8)
        for i, p in enumerate(ps):
            if round(p, precision_p) in p_done:
                print(f'Perturbation {p} has been run!')
                continue
            pool.apply_async(run, args=(max_component_subgraph, p, shortest_path_times, perturbation_times, save_path, method, rtheta_coords, net1_path, gamma, T, ), 
            callback=write_results, error_callback=print_error)
        pool.close()
        pool.join()
    else:
        for i, p in enumerate(ps):
            result = ""
            if round(p, precision_p) in p_done:
                print(f'Perturbation {p} has been run!')
                continue
            _, text_result = run(max_component_subgraph, p, shortest_path_times, perturbation_times, save_path, method, rtheta_coords, net1_path, gamma, T)
            result += text_result
            write_results((save_path, result))

def main1():
    """ reorganize main0 """
    net1_path = "./dk-lab-2020_code_hyperlink/netdata/h2v1.dat"
    my_exp = exp()
    my_exp.load_network(net_path=net1_path)
    """ prepare exp parameter """
    my_exp.precision_p = 5
    my_exp.ps = np.around(np.arange(0.1, 1, 0.001), my_exp.precision_p)
    # ps = np.concatenate((np.linspace(0.0001, 0.1, 10), np.linspace(0.5, 0.9, 2)), axis=None)
    print(f"Number of p {np.size(my_exp.ps)}")
    my_exp.shortest_path_times = 2  # 5
    my_exp.perturbation_times = 2  # 5
    # save_path = "./results/network1_overlap_results_6.24.txt"
    # save_path = "./results/network1_overlap_hgd_results_6.25.txt"
    # save_path = "./results/network1_overlap_dijkstra_results_6.25.txt"
    # save_path = f"./results/network1_overlap_dijkstra_results_{shortest_path_times}_{perturbation_times}_6.26.txt"
    my_exp.save_path = f"./results/network1_overlap_perturbate_hgd_results_{my_exp.shortest_path_times}_{my_exp.perturbation_times}_6.26.txt"
    my_exp.multiprocessing = False
    # method = "hyperbolic_geodesic_distance"
    # method = "dijkstra"
    my_exp.method = "perturbated_hyperbolic_geodesic_distance"
    if my_exp.method == "hyperbolic_geodesic_distance":
        # use given actual hyperbolic coordinates
        print(f"For method {method}, get r_theta coordinates.")
        my_exp.rtheta_coords = dict()
        with open(f"./dk-lab-2020_code_hyperlink/coords_real/h2v1.coords.dat", 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                i = int(line[0])
                r = float(line[1])
                theta = float(line[2])
                my_exp.rtheta_coords[i] = (r, theta)
    elif my_exp.method == "perturbated_hyperbolic_geodesic_distance":
        # use hyperlink to embed hyperbolic coordinates in perturbated network, given gamma, T and p
        my_exp.gamma = 2.5
        my_exp.T = 0.1
    my_exp.save_run()

def main2():
    net_path = "./dk-lab-2020_code_hyperlink/netdata/rhg_testA1_n500_alpha0.75_T0.5_d5.net"
    my_exp = exp()
    my_exp.load_network(net_path=net_path)
    """ prepare exp parameter """
    my_exp.precision_p = 5
    # my_exp.ps = np.around(np.arange(0.01, 1, 0.01), my_exp.precision_p)
    my_exp.ps = np.around(np.arange(0.05, 1, 0.05), my_exp.precision_p)  # 19
    # ps = np.concatenate((np.linspace(0.0001, 0.1, 10), np.linspace(0.5, 0.9, 2)), axis=None)
    print(f"Number of p {np.size(my_exp.ps)}")
    my_exp.shortest_path_times = 100  # 5
    my_exp.perturbation_times = 15  # 5
    # my_exp.metric = "edit_distance"
    my_exp.metric = "all"
    my_exp.multiprocessing = False
    # my_exp.method = "hyperbolic_geodesic_distance"
    # my_exp.method = "dijkstra"
    # my_exp.method = "perturbated_hyperbolic_geodesic_distance"
    # my_exp.method = "perturbated_hyperbolic_p_dijkstra"
    # my_exp.method = "all"
    my_exp.method = ["perturbated_mercator_geodesic_distance"]
    if my_exp.method == "hyperbolic_geodesic_distance":
        # use given actual hyperbolic coordinates
        print(f"For method {method}, get r_theta coordinates.")
        my_exp.rtheta_coords = dict()
        real_coords_path = f"./dk-lab-2020_code_hyperlink/coords_real/rhg_testA1_n500_alpha0.75_T0.5_d5.coords.dat"
        with open(real_coords_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                i = int(line[0])
                r = float(line[1])
                theta = float(line[2])
                my_exp.rtheta_coords[i] = (r, theta)
    elif my_exp.method[0].startswith("perturbated_hyperbolic") or my_exp.method == "all":
        # use hyperlink to embed hyperbolic coordinates in perturbated network, given gamma, T and p
        my_exp.gamma = 2.5
        my_exp.T = 0.5
        my_exp.R = 13.604
    # save_path = "./results/network1_overlap_results_6.24.txt"
    # save_path = "./results/network1_overlap_hgd_results_6.25.txt"
    # save_path = "./results/network1_overlap_dijkstra_results_6.25.txt"
    # save_path = f"./results/network1_overlap_dijkstra_results_{shortest_path_times}_{perturbation_times}_6.26.txt"
    # my_exp.save_path = f"./results/testA1_overlap_perturbate_hgd_results_{my_exp.shortest_path_times}_{my_exp.perturbation_times}_6.26.txt"
    # my_exp.save_path = f"./results/testA1_overlap_perturbate_hgd_results_{my_exp.shortest_path_times}_{my_exp.perturbation_times}_6.27.txt"
    # my_exp.save_path = f"./results/testA1_overlap_perturbate_hp_dijkstra_results_{my_exp.shortest_path_times}_{my_exp.perturbation_times}_6.28.txt"
    # my_exp.save_path = f"./results/testA1_overlap_perturbate_hgd_results_{my_exp.shortest_path_times}_{my_exp.perturbation_times}_6.28.txt"
    # my_exp.save_path = f"./results/testA1_{my_exp.metric}_dijkstra_results_{my_exp.shortest_path_times}_{my_exp.perturbation_times}_6.30.txt"
    # my_exp.save_path = f"./results/testA1_{my_exp.metric}_{my_exp.method}_results_{my_exp.shortest_path_times}_{my_exp.perturbation_times}_7.1.txt"
    # my_exp.save_path = f"./results/testA1_{my_exp.metric}_{my_exp.method}_results_{my_exp.shortest_path_times}_{my_exp.perturbation_times}_7.2.txt"
    my_exp.save_path = f"./results/testA1_{my_exp.metric}_{my_exp.method}_results_{my_exp.shortest_path_times}_{my_exp.perturbation_times}_7.14.txt"
    my_exp.save_run()

def main3():
    net_path = "./dk-lab-2020_code_hyperlink/netdata/h2v2.dat"
    my_exp = exp()
    my_exp.load_network(net_path=net_path)
    """ prepare exp parameter """
    my_exp.precision_p = 5
    # my_exp.ps = [0]
    my_exp.ps = np.around(np.arange(0.1, 1, 0.1), my_exp.precision_p)  # 9
    # ps = np.concatenate((np.linspace(0.0001, 0.1, 10), np.linspace(0.5, 0.9, 2)), axis=None)
    print(f"Number of p {np.size(my_exp.ps)}")
    my_exp.shortest_path_times = 20  # 5
    my_exp.perturbation_times = 1  # 5
    # my_exp.metric = "edit_distance"
    my_exp.metric = "all"
    my_exp.multiprocessing = False
    # my_exp.method = "hyperbolic_geodesic_distance"
    # my_exp.method = "dijkstra"
    # my_exp.method = "perturbated_hyperbolic_geodesic_distance"
    # my_exp.method = "perturbated_hyperbolic_p_dijkstra"
    # my_exp.method = "all"
    my_exp.method = ["dijkstra", "perturbated_mercator_geodesic_distance"]
    my_exp.save_path = f"./results/h2v2_{my_exp.metric}_{my_exp.method}_mercator_results_{my_exp.shortest_path_times}_{my_exp.perturbation_times}_7.14.txt"
    my_exp.save_run()


if __name__=="__main__":
    # main0()
    # main1()
    main2()
    # main3()
