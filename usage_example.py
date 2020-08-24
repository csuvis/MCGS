import os

from MCGS import MCGS
from utils.load_graph import loadGraph
from utils.save_file import saveAsSVG


def main():
    MCGS_model = MCGS()  # create an mino-centric graph sampling (MCGS) algorithm model
    dataset = 'As-733'
    G = loadGraph(dataset)  # load original graph
    Gs = MCGS_model.run_sampling(G, 0.3, alpha=1, beta=4, loss_weight=[1, 0, 0])  # run the samping algorithm
    if not os.path.exists('output'):
    	os.makedirs('output') 
    saveAsSVG(dataset, Gs)  # save the sampling graph data as SVG type file
    print('Sampling completed successfully!')
    print('Sampling Graph nodes: {}'.format(len(Gs.nodes())))
    print('Sampling Graph edges: {}'.format(len(Gs.edges())))


if __name__ == '__main__':
    main()