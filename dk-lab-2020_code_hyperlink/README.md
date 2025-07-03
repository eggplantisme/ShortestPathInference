# HyperLink: an algorithm to embed undirected unweighted networks in two-dimensional hyperbolic space

## Summary
This repository contains the code to embed networks to hyperbolic space using HyperLink algorithm as described in the paper [Link prediction with hyperbolic geometry](https://doi.org/10.1103/PhysRevResearch.2.043113). If you use this code, please cite the original paper and this repository.

## Dependencies

The embedding algorithm is contained in the _hyperlink.cpp_ file. Auxiliary functions and classes are contained in the _/auxiliary/_ directory. To compile the code, you will need:

1. g++ compiler

2. boost library

## Compilation and usage

The code is compiled using the following command:

```
g++ hyperlink.cpp auxiliary/network.cpp auxiliary/global.cpp auxiliary/mle.cpp -o hyperlink.exe -O3
```

The code then runs as follows:
```
./hyperlink.exe input_network_path gamma temperature random_seed output_coordinates_path fraction_links_removed num_layers grid_multiplier
```
Here, `input_network_path` is the input network in the edge list format, `gamma` is the power-law exponent of the degree distribution, `temperature` is the temperature parameter of the static hyperbolic random graph model ![$\mathbb{H}^2$](https://render.githubusercontent.com/render/math?math=\mathbb{H}^2&mode=inline), `random_seed` parameter provides a random seed for the random number generator, `output_coordinates_path` file stores the resulting embedding, `fraction_links_removed` is the fraction of randomly removed links in the input network (in case we use HyperLink for link prediction or other tasks with the _training_ network), `num_layers` is the number of logarithmically sized layers, and `grid_multiplier` is the coefficient controlling the size of the mesh in the angular space used to infer nodes' angles. By default, we set `num_layers` to 20, and `grid_multiplier` to 1.

**Note**: edge list should contain node pairs specified by their indices starting from 0.

## Usage example

Assume you have a network _input_network.net_ in the form of edge list that you want to embed in two-dimensional hyperbolic space. You measured that the network has the degree distribution exponent of ![$\gamma=2.5$](https://render.githubusercontent.com/render/math?math=\gamma=2.5&mode=inline), and the observed average clustering coefficient in the network gives you that the temperature parameter is ![$T=0.3$](https://render.githubusercontent.com/render/math?math=T=0.3&mode=inline). Also, assume you then want to use this embedding to perform link prediction, and you know that the training network that you provide to HyperLink has ![$p_r=0.1$](https://render.githubusercontent.com/render/math?math=p_r=0.1&mode=inline) fraction of links randomly removed from it. In this case, you should run the following code:
```
./hyperlink.exe input_network.net 2.5 0.3 random_seed output_coordinates_path 0.1 num_layers grid_multiplier
```

## HyperLink algorithm description

HyperLink is a complex embedding algorithm based on the maximum likelihood estimation (MLE) of nodes' hyperbolic coordinates. Please see the `algorithm_description.pdf` file at this repository for a detailed algorithm description and specifications.

## Edits

**Edit 1 (2021-06-02):** added an error message for the cases when finite-size parameters inference routine fails to find a solution, and updated the requirements for the edge list format. Many thanks to Filippo Radicchi and Yijiao Zhang for pointing out the issue! 

## Contact info

In case you notice any bugs or have difficulties with running the code, please feel free to contact the authors via Bitbucket or email. 

Authors of the code: Maksim Kitsak and Ivan Voitalov.
