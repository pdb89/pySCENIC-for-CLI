#!/usr/bin/env python
# coding: utf-8

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("wdir", help="working directory")
    parser.add_argument("exp_matrix", help="matrix of gene expression values exported from seurat or anndata object")
    parser.add_argument("f_tfs", help="TF database file")
    parser.add_argument("f_db_dir", help="ranking database (feather files) directory")
    parser.add_argument("f_motif_path", help="motif database file from tf2motif")
    parser.add_argument("outdir", help="output directory")
    args = parser.parse_args()



    # import dependencies
    import os
    import anndata
    import pickle
    import gc
    import glob
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import loompy as lp
    import matplotlib.pyplot as plt


    from MulticoreTSNE import MulticoreTSNE as TSNE
    from dask.diagnostics import ProgressBar

    from arboreto.utils import load_tf_names
    from arboreto.algo import grnboost2

    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
    from pyscenic.utils import modules_from_adjacencies, load_motifs
    from pyscenic.prune import prune2df, df2regulons


    ## ISSUE with Numpy version (np.object - deprecated after v1.20, currently v2.1.3 which is req'd for other pkgs)
    ## TEMP FIX: manually changed "np.object" to "object" in ~/tensorflow_venv/lib/python3.11/site-packages/pyscenic/transform.py
    ## (backup of original file in ~/)
    from pyscenic.prune import prune2df, df2regulons

    from pyscenic.aucell import aucell

    import seaborn as sns

    # ## SCENIC steps - starting with expression matrix (export from Seurat or adata object, i.e.)
    # ### STEP 1: Gene regulatory network inference, and generation of co-expression modules
    # #### Phase Ia: GRN inference using the GRNBoost2 algorithm
    # 
    # For this step the CLI version of SCENIC is used. This step can be deployed on an High Performance Computing system. We use the counts matrix (without log transformation or further processing) from the loom file we wrote earlier.
    # _Output:_ List of adjacencies between a TF and its targets stored in `ADJACENCIES_FNAME`.


    # set variables for file paths to read from and write to:

    # set a working directory
    wdir = args.wdir  #"/storage/home/bou01pav/Projects/Saul_lab/pySCENIC/pbmc3k_TESTING/"
    os.chdir( wdir )

    # set directory for output files
    outdir = args.outdir

    # path to unfiltered loom file (this will be created in the optional steps below)
    f_loom_path_unfilt = outdir + "/expr_pbmc3k_unfiltered.loom" # test dataset, n=500 cells

    # # path to loom file with basic filtering applied (this will be created in the "initial filtering" step below). Optional.
    f_loom_path_scenic = outdir + "/expr_pbmc3k_filtered_scenic.loom"

    ## # path to anndata object, which will be updated to store Scanpy results as they are generated below
    #f_anndata_path = wdir + "pbmc3k.h5ad"

    # path to pyscenic output
    f_pyscenic_output = outdir + "/pyscenic_output.loom"

    # loom output, generated from a combination of Scanpy and pySCENIC results:
    f_final_loom = outdir + '/expr_pbmc3k_scenic_integrated-output.loom'

    # housekeeping
    sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
    #sc.logging.print_versions()
    sc.set_figure_params(dpi=150, fontsize=10, dpi_save=600)


    # #### Load in the expression matrix (note the delimiter: tsv or csv)
    f_exprMat = args.exp_matrix #"/expr_mini_pbmc3k.tsv"
    ex_matrix = pd.read_csv(f_exprMat, sep = "\t", header=0, index_col=0).T  # TRANSPOSED expression matrix (rows = cells, cols = genes)

    # #### OR load in anndata file (exported/converted from saved Seurat)
    #adata = anndata.io.read_h5ad(f_anndata_path)

    ## restrict expression matrix to ONLY genes present in the .feather file
    #ranking_feather = pd.read_feather("/storage/home/bou01pav/Reference/SCENIC/cisTarget_databases/hg19_refseqr45_v9_genebased/hg19-500bp-upstream-10species.mc9nr.genes_vs_motifs.rankings.feather")
    #overlap_values = ex_matrix.index[pd.Series(ex_matrix.index).isin(ranking_feather.columns)].unique()
    #ex_matrix = ex_matrix.loc[overlap_values, :].T


    # # create basic row and column attributes for the loom file:
    # row_attrs = {
    #     "Gene": np.array(adata.var_names) ,
    # }
    # col_attrs = {
    #     "CellID": np.array(adata.obs_names) ,
    #     "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    #     "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
    # }
    # lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)


    # #### Load in the TF list file

    # transcription factors list
    f_tfs = args.f_tfs  #"/storage/home/bou01pav/Projects/Saul_lab/pySCENIC/pbmc3k_TESTING/ref/hs_hgnc_tfs.txt" 
    tf_names = load_tf_names( f_tfs )


    # #### List and load in the pre-downloaded ranking database (.feather) files and motif database file

    # ranking databases
    f_db_glob = args.f_db_dir + "/*feather"
    f_db_names = ' '.join( glob.glob(f_db_glob) )

    # motif databases
    f_motif_path = args.f_motif_path  #"/storage/home/bou01pav/Reference/SCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

    db_fnames = glob.glob(f_db_glob)
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
    dbs


    # #### Get list of adjacencies between TFs and targets
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
    gc.collect()
    #adjacencies.head(n = 10)

    # Write adjacencies matrix to file
    adj_out_path = args.outdir + "/adj.tsv"
    adjacencies.to_csv(adj_out_path, sep="\t")


    # ### STEP 2-3: Regulon prediction aka cisTarget from CLI
    # 
    # _Output:_ List of adjacencies between a TF and its targets stored in `motifs_fname`.

    # locations for ranking databases, and motif annotations:

    # Here, we use the `--mask_dropouts` option, which affects how the correlation between TF and target genes is calculated during module creation. It is important to note that prior to pySCENIC v0.9.18, the default behavior was to mask dropouts, while in v0.9.18 and later, the correlation is performed using the entire set of cells (including those with zero expression). When using the `modules_from_adjacencies` function directly in python instead of via the command line, the `rho_mask_dropouts` option can be used to control this.


    # create TF modules (with dropout masking)
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix, rho_mask_dropouts=True))
    gc.collect()

    # paths for enriched motif and regulon output files
    motifs_fname = args.outdir + "/motifs.csv"  #"/storage/home/bou01pav/Projects/Saul_lab/pySCENIC/pbmc3k_TESTING/motifs.csv"
    regulons_name = args.outdir + "/regulons.p"  #"/storage/home/bou01pav/Projects/Saul_lab/pySCENIC/pbmc3k_TESTING/regulons.p"


    # Calculate a list of enriched motifs and the corresponding target genes for all modules.
    with ProgressBar():
        df = prune2df(dbs, modules, f_motif_path, filter_for_annotation=True,
                      client_or_address = 'custom_multiprocessing',
                      num_workers = len(db_fnames))
    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)
    # Save the enriched motifs and the discovered regulons to disk.
    df.to_csv(motifs_fname)
    with open(regulons_name, "wb") as f:
        pickle.dump(regulons, f)
    gc.collect()


    # Reload the enriched motifs and regulons from file
    df = load_motifs(motifs_fname)
    with open(regulons_name, "rb") as f:
        regulons = pickle.load(f)
    gc.collect()


    # ### STEP 4: Cellular enrichment (aka AUCell)
    # It is important to check that most cells have a substantial fraction of expressed/detected genes in the calculation of the AUC. 

    # #### Calculate signature enrichment scores on each cell
    auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
    g = sns.clustermap(auc_mtx, figsize=(8,8))
    g.fig.savefig("clustermap.png")
    gc.collect()

    # save the AUC matrix
    auc_mtx_outpath = args.outdir + "/auc_mtx.csv"
    auc_mtx.to_csv(auc_mtx_outpath)
    gc.collect()

    nGenesDetectedPerCell = np.sum(ex_matrix>0, axis=1)
    percentiles = nGenesDetectedPerCell.quantile([.01, .05, .10, .50, 1])
    print(percentiles)


    # The following histogram gives an idea of the distribution and allows selection of an appropriate threshold.
    # In this plot, a few thresholds are highlighted, with the number of genes selected shown in red text and the corresponding percentile in parentheses).
    # See [the relevant section in the R tutorial](https://bioconductor.org/packages/devel/bioc/vignettes/AUCell/inst/doc/AUCell.html#build-gene-expression-rankings-for-each-cell) for more information.

    fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=150)
    sns.distplot(nGenesDetectedPerCell, norm_hist=False, kde=False, bins='fd')
    for i,x in enumerate(percentiles):
        fig.gca().axvline(x=x, ymin=0,ymax=1, color='red')
        ax.text(x=x, y=ax.get_ylim()[1], s=f'{int(x)} ({percentiles.index.values[i]*100}%)', color='red', rotation=30, size='x-small',rotation_mode='anchor' )
    ax.set_xlabel('# of genes')
    ax.set_ylabel('# of cells')
    fig.tight_layout() 
    fig_outpath = args.outdir + "/histogram.png"
    fig.savefig(fig_outpath)
    gc.collect()


    # ### Visualization of SCENIC's AUC matrix

    # First, load the relevant data from the loom we just created

    import json
    import zlib
    import base64
    import umap

    ## collect SCENIC AUCell output
    #lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
    #auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
    #lf.close()


    # UMAP
    runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
    dr_umap = runUmap( auc_mtx )
    umap_outpath = args.outdir + "/scenic_umap.txt"
    pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv(umap_outpath, sep='\t')
    gc.collect()

    # tSNE
    tsne = TSNE( n_jobs=20 )
    dr_tsne = tsne.fit_transform( auc_mtx )
    tsne_outpath = args.outdir + "/scenic_tsne.txt"
    pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv(tsne_outpath, sep='\t')
    gc.collect()


    # ## Integrate the output
    # 
    # Here, we combine the results from SCENIC and the Scanpy analysis into a SCope-compatible loom file

    # In[ ]:


    ## scenic output
    #lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
    #meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
    #exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID)
    #auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
    #regulons = lf.ra.Regulons
    #dr_umap = pd.read_csv( 'scenic_umap.txt', sep='\t', header=0, index_col=0 )
    #dr_tsne = pd.read_csv( 'scenic_tsne.txt', sep='\t', header=0, index_col=0 )
    ###


    ## Fix regulon objects to display properly in SCope:
    #auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
    #regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
    ## regulon thresholds
    #rt = meta['regulonThresholds']
    #for i,x in enumerate(rt):
    #    tmp = x.get('regulon').replace("(","_(")
    #    x.update( {'regulon': tmp} )


    # Concatenate embeddings (tSNE, UMAP, etc.)

    # In[ ]:


    #tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])

    #Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
    #Embeddings_X = pd.concat( [
    #        pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[0] ,
    #        pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0] ,
    #        dr_tsne['X'] ,
    #        dr_umap['X']
    #    ], sort=False, axis=1, join='outer' )
    #Embeddings_X.columns = ['1','2','3','4']

    #Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
    #Embeddings_Y = pd.concat( [
    #        pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[1] ,
    #        pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1] ,
    #        dr_tsne['Y'] ,
    #        dr_umap['Y']
    #    ], sort=False, axis=1, join='outer' )
    #Embeddings_Y.columns = ['1','2','3','4']


    # Metadata:

    # In[ ]:


    ### metadata
    # metaJson = {}

    # metaJson['embeddings'] = [
    #     {
    #         "id": -1,
    #         "name": f"Scanpy t-SNE (highly variable genes)"
    #     },
    #     {
    #         "id": 1,
    #         "name": f"Scanpy UMAP  (highly variable genes)"
    #     },
    #     {
    #         "id": 2,
    #         "name": "Scanpy PC1/PC2"
    #     },
    #     {
    #         "id": 3,
    #         "name": "SCENIC AUC t-SNE"
    #     },
    #     {
    #         "id": 4,
    #         "name": "SCENIC AUC UMAP"
    #     },
    # ]

    # metaJson["clusterings"] = [{
    #             "id": 0,
    #             "group": "Scanpy",
    #             "name": "Scanpy louvain default resolution",
    #             "clusters": [],
    #         }]

    # metaJson["metrics"] = [
    #         {
    #             "name": "nUMI"
    #         }, {
    #             "name": "nGene"
    #         }, {
    #             "name": "Percent_mito"
    #         }
    # ]

    # metaJson["annotations"] = [
    #     {
    #         "name": "Louvain_clusters_Scanpy",
    #         "values": list(set( adata.obs['louvain'].astype(np.str) ))
    #     },
    #     #{
    #     #    "name": "Genotype",
    #     #    "values": list(set(adata.obs['Genotype'].values))
    #     #},
    #     #{
    #     #    "name": "Timepoint",
    #     #    "values": list(set(adata.obs['Timepoint'].values))
    #     #},
    #     #{
    #     #    "name": "Sample",
    #     #    "values": list(set(adata.obs['Sample'].values))
    #     #}
    # ]

    # # SCENIC regulon thresholds:
    # metaJson["regulonThresholds"] = rt

    # for i in range(max(set([int(x) for x in adata.obs['louvain']])) + 1):
    #     clustDict = {}
    #     clustDict['id'] = i
    #     clustDict['description'] = f'Unannotated Cluster {i + 1}'
    #     metaJson['clusterings'][0]['clusters'].append(clustDict)
        
    # clusterings = pd.DataFrame()
    # clusterings["0"] = adata.obs['louvain'].values.astype(np.int64)


    # Assemble loom file row and column attributes

    # In[ ]:


    # def dfToNamedMatrix(df):
    #     arr_ip = [tuple(i) for i in df.values]
    #     dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    #     arr = np.array(arr_ip, dtype=dtyp)
    #     return arr


    # In[ ]:


    # col_attrs = {
    #     "CellID": np.array(adata.obs.index),
    #     "nUMI": np.array(adata.obs['n_counts'].values),
    #     "nGene": np.array(adata.obs['n_genes'].values),
    #     "Louvain_clusters_Scanpy": np.array( adata.obs['louvain'].values ),
    #     #"Genotype": np.array(adata.obs['Genotype'].values),
    #     #"Timepoint": np.array(adata.obs['Timepoint'].values),
    #     #"Sample": np.array(adata.obs['Sample'].values),
    #     "Percent_mito": np.array(adata.obs['percent_mito'].values),
    #     "Embedding": dfToNamedMatrix(tsneDF),
    #     "Embeddings_X": dfToNamedMatrix(Embeddings_X),
    #     "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
    #     "RegulonsAUC": dfToNamedMatrix(auc_mtx),
    #     "Clusterings": dfToNamedMatrix(clusterings),
    #     "ClusterID": np.array(adata.obs['louvain'].values)
    # }

    # row_attrs = {
    #     "Gene": lf.ra.Gene,
    #     "Regulons": regulons,
    # }

    # attrs = {
    #     "title": "sampleTitle",
    #     "MetaData": json.dumps(metaJson),
    #     "Genome": 'hg38',
    #     "SCopeTreeL1": "",
    #     "SCopeTreeL2": "",
    #     "SCopeTreeL3": ""
    # }

    # # compress the metadata field:
    # attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')


    # Create a new loom file, copying the expression matrix from the open loom connection:

    # In[ ]:


    # lp.create(
    #     filename = f_final_loom ,
    #     layers=lf[:,:],
    #     row_attrs=row_attrs, 
    #     col_attrs=col_attrs, 
    #     file_attrs=attrs
    # )
    # lf.close() # close original pyscenic loom file


    # This loom file can now be imported into [SCope](http://scope.aertslab.org/).

if __name__ == "__main__":
    main()
