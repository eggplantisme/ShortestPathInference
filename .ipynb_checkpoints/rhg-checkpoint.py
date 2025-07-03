import scipy.stats as st
import numpy as np
import networkx as nx
from tqdm import tqdm


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
        for i in range(self.n):
            for j in range(i + 1, self.n):
                delta_theta = np.pi-np.abs(np.pi-np.abs(thetas[i]-thetas[j]))
                d = np.arccosh(np.cosh(rs[i])*np.cosh(rs[j])-np.sinh(rs[i])*np.sinh(rs[j])*np.cos(delta_theta))
                X[i, j] = d
                X[j, i] = d
                p = 1/(np.exp((d-self.R)/self.T)+1)
                A[i, j] = A[j, i] = np.random.binomial(1, p)
        return A
    
    @staticmethod
    def save(A, path):
        G = nx.from_numpy_array(A)
        with open(path, 'w') as f:
            for line in nx.generate_edgelist(G, data=False):
                f.write(line+'\n')

    @staticmethod
    def load(path):
        edges = []
        with open(path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip().split(' ')
                edges.append((int(line[0]), int(line[1])))
        G = nx.from_edgelist(edges, create_using=nx.Graph())
        return G



class r_pdf(st.rv_continuous):
    def __init__(self, a, b, alpha, R, name='r_pdf'):
        super().__init__(a=a, b=b, name=name)
        self.alpha = alpha
        self.R = R

    def _pdf(self,x):
        return np.sinh(self.alpha*x)/(np.cosh(self.alpha*self.R)-1)