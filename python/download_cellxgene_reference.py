#!/usr/bin/env python3
"""
download_cellxgene_reference.py
===============================================================================
Download reference datasets from CELLxGENE Census for ACP analysis.

Downloads skin keratinocytes and/or oral mucosa epithelium for comparison
with ACP epithelial subtypes.

Usage:
    python scripts/download_cellxgene_reference.py
    python scripts/download_cellxgene_reference.py --tissue skin --max_cells 50000
    python scripts/download_cellxgene_reference.py --tissue oral
    python scripts/download_cellxgene_reference.py --tissue all

Output:
    data/reference/cellxgene_skin_keratinocytes.h5ad
    data/reference/cellxgene_oral_epithelium.h5ad
===============================================================================
"""

import argparse
import os
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="Download reference data from CELLxGENE Census"
    )
    parser.add_argument(
        "--tissue", "-t",
        choices=["skin", "oral", "all"],
        default="skin",
        help="Tissue to download (default: skin)"
    )
    parser.add_argument(
        "--output_dir", "-o",
        default="data/reference",
        help="Output directory (default: data/reference)"
    )
    parser.add_argument(
        "--max_cells", "-m",
        type=int,
        default=100000,
        help="Maximum cells to download per tissue (default: 100000)"
    )
    parser.add_argument(
        "--healthy_only",
        action="store_true",
        default=True,
        help="Only download healthy/normal samples (default: True)"
    )
    parser.add_argument(
        "--census_version",
        default="stable",
        help="CELLxGENE Census version (default: stable)"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("\n" + "=" * 64)
    print("CELLxGENE Census Reference Data Download")
    print("=" * 64)
    print(f"Tissue: {args.tissue}")
    print(f"Max cells: {args.max_cells:,}")
    print(f"Healthy only: {args.healthy_only}")
    print(f"Output: {args.output_dir}")
    print(f"Census version: {args.census_version}")
    print("=" * 64 + "\n")
    
    # Import after arg parsing for faster --help
    try:
        import cellxgene_census
        import scanpy as sc
        import pandas as pd
        import numpy as np
        print(f"cellxgene_census: {cellxgene_census.__version__}")
        print(f"scanpy: {sc.__version__}")
    except ImportError as e:
        print(f"Error: {e}")
        print("\nPlease install required packages:")
        print("  pip install cellxgene-census scanpy anndata")
        sys.exit(1)
    
    # Define queries
    tissue_queries = {
        "skin": {
            "name": "Skin Keratinocytes",
            "filter": (
                "tissue_general == 'skin of body' and "
                "is_primary_data == True and "
                "(cell_type == 'keratinocyte' or "
                "cell_type == 'basal cell of epidermis' or "
                "cell_type == 'spinous cell of epidermis')"
            ),
            "output_file": "cellxgene_skin_keratinocytes.h5ad"
        },
        "oral": {
            "name": "Oral Mucosa Epithelium",
            "filter": (
                "(tissue_general == 'tongue' or "
                "tissue_general == 'mouth' or "
                "tissue == 'gingiva' or "
                "tissue == 'oral mucosa' or "
                "tissue == 'buccal mucosa') and "
                "is_primary_data == True and "
                "(cell_type == 'keratinocyte' or "
                "cell_type == 'epithelial cell' or "
                "cell_type == 'oral epithelial cell' or "
                "cell_type == 'basal cell')"
            ),
            "output_file": "cellxgene_oral_epithelium.h5ad"
        }
    }
    
    # Add healthy filter
    if args.healthy_only:
        for tissue in tissue_queries:
            tissue_queries[tissue]["filter"] += " and disease == 'normal'"
    
    # Determine what to download
    tissues = list(tissue_queries.keys()) if args.tissue == "all" else [args.tissue]
    
    # Open Census connection
    print("Opening CELLxGENE Census...")
    census = cellxgene_census.open_soma(census_version=args.census_version)
    
    try:
        for tissue in tissues:
            query = tissue_queries[tissue]
            output_path = os.path.join(args.output_dir, query["output_file"])
            
            print(f"\n--- Downloading: {query['name']} ---")
            
            # Check if exists
            if os.path.exists(output_path):
                print(f"File already exists: {output_path}")
                print("Delete to re-download.")
                continue
            
            # Count cells
            print("Querying cell count...")
            try:
                obs_df = cellxgene_census.get_obs(
                    census,
                    organism="Homo sapiens",
                    value_filter=query["filter"],
                    column_names=["cell_type", "tissue", "disease"]
                )
                total_cells = len(obs_df)
                print(f"Found {total_cells:,} cells matching query")
                
                if total_cells == 0:
                    print("No cells found. Trying broader query...")
                    # Try without cell type filter - just get all skin cells
                    broader_filter = "tissue_general == 'skin of body' and is_primary_data == True"
                    if args.healthy_only:
                        broader_filter += " and disease == 'normal'"
                    obs_df = cellxgene_census.get_obs(
                        census,
                        organism="Homo sapiens",
                        value_filter=broader_filter,
                        column_names=["cell_type", "tissue", "disease"]
                    )
                    print(f"Broader query found {len(obs_df):,} cells")
                    print("Cell types available:")
                    print(obs_df["cell_type"].value_counts().head(20))
                    continue
                    
            except Exception as e:
                print(f"Query failed: {e}")
                continue
            
            # Subsample if needed
            obs_coords = None
            if total_cells > args.max_cells:
                print(f"Subsampling to {args.max_cells:,} cells...")
                sample_idx = np.random.choice(total_cells, args.max_cells, replace=False)
                obs_coords = sample_idx.tolist()
            
            # Download
            print("Downloading gene expression data...")
            print("(This may take several minutes...)")
            
            adata = cellxgene_census.get_anndata(
                census,
                organism="Homo sapiens",
                obs_value_filter=query["filter"],
                obs_column_names=[
                    "cell_type",
                    "cell_type_ontology_term_id", 
                    "tissue",
                    "tissue_general",
                    "disease",
                    "assay",
                    "donor_id",
                    "dataset_id",
                    "sex",
                    "development_stage"
                ],
                obs_coords=obs_coords
            )
            
            print(f"Downloaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
            
            # Save
            print(f"Saving to: {output_path}")
            adata.write_h5ad(output_path)
            
            # Print summary
            print("\nCell type distribution:")
            print(adata.obs["cell_type"].value_counts().head(10))
            
            print("\nTissue distribution:")
            print(adata.obs["tissue"].value_counts().head(10))
            
            # File size
            size_mb = os.path.getsize(output_path) / 1e6
            print(f"\nFile size: {size_mb:.1f} MB")
            
    finally:
        census.close()
    
    print("\n" + "=" * 64)
    print("Download Complete")
    print("=" * 64)
    print(f"\nFiles saved to: {args.output_dir}/")
    print("\nTo load in R:")
    print("  library(anndata)")
    print(f"  adata <- read_h5ad('{args.output_dir}/cellxgene_skin_keratinocytes.h5ad')")
    print("\nTo load in Python:")
    print("  import scanpy as sc")
    print(f"  adata = sc.read_h5ad('{args.output_dir}/cellxgene_skin_keratinocytes.h5ad')")
    print()


if __name__ == "__main__":
    main()
