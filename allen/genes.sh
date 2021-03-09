#!/bin/sh

wget "http://api.brain-map.org/api/v2/data/query.csv?criteria=model::Gene,rma::criteria,products[abbreviation\$eq'DevMouse'],rma::options,[tabular\$eq'genes.id','genes.acronym+as+gene_symbol','genes.name+as+gene_name','genes.entrez_id+as+entrez_gene_id','genes.homologene_id+as+homologene_group_id'],[order\$eq'genes.acronym']&num_rows=all&start_row=0"
