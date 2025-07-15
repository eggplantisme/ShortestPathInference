import scipy.stats as st
import numpy as np
import networkx as nx
from tqdm import tqdm
import re
from scipy.sparse import csr_matrix, lil_array
import os
# import mercator


class rhg:
    def __init__(self, n, R, T, alpha):
        self.n = n
        self.R = R
        self.T = T
        self.alpha = alpha
        pass

    def generate(self):
        """
        Generate a random hyperbolic graph based on the parameters.
        """
        # sample radii and angles
        rs = []
        thetas = []
        for i in range(self.n):
            rs.append(r_pdf(a=0, b=self.R, alpha=self.alpha, R=self.R).rvs())
            thetas.append(np.random.uniform(0, 2 * np.pi))
        # calculate the distance matrix
        X = np.zeros((self.n, self.n))
        A = np.zeros((self.n, self.n))
        for i in tqdm(range(self.n)):
            for j in range(i + 1, self.n):
                delta_theta = np.pi-np.abs(np.pi-np.abs(thetas[i]-thetas[j]))
                d = np.arccosh(np.cosh(rs[i])*np.cosh(rs[j])-np.sinh(rs[i])*np.sinh(rs[j])*np.cos(delta_theta))
                X[i, j] = d
                X[j, i] = d
                p = 1/(np.exp((d-self.R)/self.T)+1)
                A[i, j] = A[j, i] = np.random.binomial(1, p)
        return A, rs, thetas
    
    @staticmethod
    def hyperlink_embedding(network_path, gamma, T, ouput_path, p=0, node_list=None, remove_output_file=False):
        command = f"./dk-lab-2020_code_hyperlink/hyperlink_jiaze.exe {network_path} {gamma} {T} 1 {ouput_path} {p} 20 1"
        print(f"command: {command}")
        os.system(command)
        rtheta_coords = dict()
        with open(ouput_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                i = int(line[0])
                r = float(line[2])
                theta = float(line[3])
                if node_list is None:
                    rtheta_coords[i] = (r, theta)
                else:
                    rtheta_coords[node_list[i]] = (r, theta)
        if remove_output_file:
            os.remove(ouput_path)
        return rtheta_coords
    
    @staticmethod
    def mercator_embedding(network_path, output_path, remove_output_file=False):
        mercator.embed(network_path, output_name=output_path)
        r_theta_coords = dict()
        with open(output_path+".inf_coord", 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("#"):
                    continue
                line = line.strip().split()
                i = int(line[0])
                r = float(line[3])
                theta = float(line[2])
                r_theta_coords[i] = (r, theta)
        if remove_output_file:
            os.remove(output_path+".inf_coord")
        return r_theta_coords

    @staticmethod
    def distance(rtheta0, rtheta1):
        delta_theta = np.pi - np.abs(np.pi - np.abs(rtheta0[1]-rtheta1[1]))
        dis = np.arccosh(np.cosh(rtheta0[0])*np.cosh(rtheta1[0])-np.sinh(rtheta0[0])*np.sinh(rtheta1[0])*np.cos(delta_theta))
        return dis
    
    @staticmethod
    def save(A, path):
        G = nx.from_numpy_array(A)
        with open(path, 'w') as f:
            for line in nx.generate_edgelist(G, data=False):
                f.write(line+'\n')
    
    @staticmethod
    def save_byedges(edges, path):
        with open(path, 'w') as f:
            for edge in edges:
                f.write(f"{edge[0]} {edge[1]}"+'\n')

    @staticmethod
    def save_coords(rs, thetas, path):
        with open(path, 'w') as f:
            for i in range(len(rs)):
                f.write(f"{i}\t{rs[i]}\t{thetas[i]}\n")


    @staticmethod
    def load(path):
        edges = []
        with open(path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = re.split(' |\t', line)
                # line = line.strip().split(' ')
                edges.append((int(line[0]), int(line[1])))
        G = nx.from_edgelist(edges, create_using=nx.Graph())
        return G
    
    @staticmethod
    def load_hyperlink_coords(path):
        rtheta_coords = dict()
        with open(path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split('\t')
                i = int(line[0])
                r = float(line[2])
                theta = float(line[3])
                rtheta_coords[i] = (r, theta)
        return rtheta_coords
    
    @staticmethod
    def random_perturbation_byA(A, p=0.1):
        """
        Randomly perturb (delete) the adjacency matrix A with probability p.
        """
        n = A.shape[0]
        perturbed_A = lil_array(A.copy())
        perturbed_A[np.tril_indices(n)] = 0
        mask = (perturbed_A.toarray() == 1) & (np.random.rand(n, n) < p)
        perturbed_A[mask] = 0
        perturbed_A = perturbed_A.tocsr()
        perturbed_A = perturbed_A + perturbed_A.T
        return perturbed_A
    
    @staticmethod
    def random_perturbation_byedges(edges, p=0.1):
        """
        Randomly perturb (delete) edges with probability p.
        """
        perturbed_edges = []
        for edge in edges:
            if np.random.rand() > p:
                perturbed_edges.append(edge)
        return perturbed_edges



class r_pdf(st.rv_continuous):
    def __init__(self, a, b, alpha, R, name='r_pdf'):
        super().__init__(a=a, b=b, name=name)
        self.alpha = alpha
        self.R = R

    def _pdf(self,x):
        return np.sinh(self.alpha*x)/(np.cosh(self.alpha*self.R)-1)